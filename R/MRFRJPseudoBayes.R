#' @title Reversible-Jump sampler for MRF interaction structure and parameters
#'
#' @description
#' An R6 class extending [MRFPseudoBayes] implementing a Reversible-Jump
#' Metropolis-Hastings algorithm for simultaneous Bayesian inference of the
#' MRF parameter array (`theta`) and the interaction structure (`mrfi`).
#'
#' @details
#' The algorithm explores a union of parameter spaces indexed by which subsets
#' of positions from a user-supplied **maximal** `mrfi` are active. A position
#' is *active* when its corresponding `theta` entries are non-zero; the current
#' active set is tracked via an `included` logical vector of length equal to
#' the number of positions in the maximal `mrfi`.
#'
#' At each iteration one of five moves is proposed:
#' \describe{
#'   \item{`within`}{Random-walk perturbation of the active `theta` entries.}
#'   \item{`birth`}{Add a randomly chosen inactive position (draw its `theta`
#'     from a `N(0, sdbirth^2)` proposal, with a conservative redistribution
#'     of that value back to currently active positions).}
#'   \item{`death`}{Remove a randomly chosen active position (set its `theta`
#'     to zero, redistributing the removed mass to the remaining active
#'     positions).}
#'   \item{`swap`}{Move the `theta` values from one active position to a
#'     randomly chosen inactive one (model dimension unchanged).}
#'   \item{`jump`}{Toggle a randomly chosen position: if active, zero it out;
#'     if inactive, add it with a fresh `N(0, sdbirth^2)` draw.}
#' }
#'
#' The log pseudo-posterior target is:
#' \deqn{\log \tilde\pi(\theta, M | z) =
#'   \log PL(z|\theta) + \sum_{k \in M} \log\phi(\theta_k; 0, \sigma_{\rm prior})
#'   - \lambda |M|}
#' where \eqn{M} is the active set, \eqn{|M|} is its cardinality, and
#' \eqn{\lambda = \log(nm)} is a BIC-like penalty.
#'
#' @importFrom mrf2d expand_array smr_stat smr_array fit_pl pl_mrf2d
#'   vec_description mrfi_to_string
#' @importFrom tidyr pivot_longer
#' @importFrom tibble as_tibble
#' @importFrom ggplot2 ggplot aes geom_line geom_hline theme_bw facet_wrap
#' @importFrom dplyr group_by summarize mutate ungroup filter
#' @importFrom stats dnorm rnorm runif quantile sd rgamma
#' @importFrom glue glue
#'
#' @export
MRFRJPseudoBayes <- R6::R6Class(
  classname = "MRFRJPseudoBayes",
  inherit = MRFPseudoBayes,

  private = list(
    .included       = NULL,   # logical vector, length npos
    .npos           = NULL,   # integer: number of positions in maximal mrfi
    .dim_per_group  = NULL,   # integer: theta entries per position
    .logpenalty     = NULL,   # numeric: BIC-like penalty per active position
    .sdbirth        = NULL,   # numeric: sd of birth proposal
    .kernel_probs   = NULL,   # numeric[5]: (within,swap,death,birth,jump) weights
    .included_chain = NULL,   # logical matrix: nsamples x npos
    .t_chain        = NULL,   # integer vector: iteration index per row (negative = warmup)

    # ------------------------------------------------------------------
    # Next block of iteration indices: `warmup` steps get negative indices
    # continuing from the most negative index used so far, and `nsamples`
    # steps get positive indices continuing from the largest used so far.
    # ------------------------------------------------------------------
    next_t = function(warmup, nsamples) {
      t_chain  <- private$.t_chain
      min_t    <- if (length(t_chain) == 0L) 0L else min(0L, min(t_chain))
      max_t    <- if (length(t_chain) == 0L) 0L else max(0L, max(t_chain))
      warmup_t <- if (warmup   > 0L) seq(min_t - warmup, min_t - 1L) else integer(0)
      main_t   <- if (nsamples > 0L) seq(max_t + 1L, max_t + nsamples) else integer(0)
      c(warmup_t, main_t)
    },

    # ------------------------------------------------------------------
    # Move: within — perturb active theta entries only
    # ------------------------------------------------------------------
    step_within = function(theta, lpl) {
      dpg  <- private$.dim_per_group
      mask <- as.numeric(rep(private$.included, each = dpg))

      proposed     <- theta + rnorm(private$.fdim, 0, private$.sdkernel) * mask
      proposed_lpl <- private$log_pl(proposed)

      logA <- (proposed_lpl - lpl) +
        sum(dnorm(proposed[proposed != 0], sd = private$.sdprior, log = TRUE)) -
        sum(dnorm(theta[theta != 0],       sd = private$.sdprior, log = TRUE))

      if (log(runif(1)) < logA)
        list(theta = proposed, lpl = proposed_lpl)
      else
        list(theta = theta, lpl = lpl)
    },

    # ------------------------------------------------------------------
    # Move: birth — add one inactive position
    # ------------------------------------------------------------------
    step_birth = function(theta, lpl) {
      dpg      <- private$.dim_per_group
      included <- private$.included
      n_incl   <- sum(included)
      npos     <- private$.npos

      if (n_incl >= npos) return(list(theta = theta, lpl = lpl))

      pos_to_add   <- sample(which(!included), 1L)
      theta_to_add <- rnorm(dpg, 0, private$.sdbirth)

      proposed <- theta
      idx <- ((pos_to_add - 1L) * dpg + 1L):(pos_to_add * dpg)
      proposed[idx] <- theta_to_add

      # Redistribute theta_to_add back to currently active positions
      if (n_incl > 0L) {
        wts  <- rgamma(n_incl, shape = 0.1)
        wts  <- wts / sum(wts)
        mult <- as.numeric(rep(included, each = dpg))
        mult[as.logical(mult)] <- rep(wts, each = dpg)
        proposed <- proposed - theta_to_add * mult
      }

      proposed_lpl <- private$log_pl(proposed)

      logA <- (proposed_lpl - lpl) +
        sum(dnorm(proposed[proposed != 0], sd = private$.sdprior, log = TRUE)) -
        sum(dnorm(theta[theta != 0],       sd = private$.sdprior, log = TRUE)) -
        private$.logpenalty -
        sum(dnorm(theta_to_add, sd = private$.sdbirth, log = TRUE)) +
        log(private$.kernel_probs[3L] / (n_incl + 1L)) -
        log(private$.kernel_probs[4L] / (npos - n_incl))

      if (log(runif(1)) < logA) {
        private$.included[pos_to_add] <- TRUE
        list(theta = proposed, lpl = proposed_lpl)
      } else {
        list(theta = theta, lpl = lpl)
      }
    },

    # ------------------------------------------------------------------
    # Move: death — remove one active position
    # ------------------------------------------------------------------
    step_death = function(theta, lpl) {
      dpg      <- private$.dim_per_group
      included <- private$.included
      n_incl   <- sum(included)
      npos     <- private$.npos

      if (n_incl == 0L) return(list(theta = theta, lpl = lpl))

      pos_to_remove  <- sample(which(included), 1L)
      idx_rem        <- ((pos_to_remove - 1L) * dpg + 1L):(pos_to_remove * dpg)
      theta_to_remove <- theta[idx_rem]

      remaining <- included
      remaining[pos_to_remove] <- FALSE

      proposed      <- theta
      proposed[idx_rem] <- 0

      # Redistribute removed theta mass to remaining active positions
      if (sum(remaining) > 0L) {
        wts  <- rgamma(sum(remaining), shape = 0.1)
        wts  <- wts / sum(wts)
        mult <- as.numeric(rep(remaining, each = dpg))
        mult[as.logical(mult)] <- rep(wts, each = dpg)
        proposed <- proposed + theta_to_remove * mult
      }

      proposed_lpl <- private$log_pl(proposed)

      logA <- (proposed_lpl - lpl) +
        sum(dnorm(proposed[proposed != 0], sd = private$.sdprior, log = TRUE)) -
        sum(dnorm(theta[theta != 0],       sd = private$.sdprior, log = TRUE)) +
        private$.logpenalty +
        sum(dnorm(theta_to_remove, sd = private$.sdbirth, log = TRUE)) +
        log(private$.kernel_probs[4L] / (npos - n_incl + 1L)) -
        log(private$.kernel_probs[3L] / n_incl)

      if (log(runif(1)) < logA) {
        private$.included[pos_to_remove] <- FALSE
        list(theta = proposed, lpl = proposed_lpl)
      } else {
        list(theta = theta, lpl = lpl)
      }
    },

    # ------------------------------------------------------------------
    # Move: swap — move theta values from one active to one inactive position
    # ------------------------------------------------------------------
    step_swap = function(theta, lpl) {
      dpg      <- private$.dim_per_group
      included <- private$.included
      n_incl   <- sum(included)
      npos     <- private$.npos

      # Need at least 2 included so one survives AND at least 1 excluded
      if (n_incl <= 1L || n_incl >= npos) return(list(theta = theta, lpl = lpl))

      to_remove <- sample(which( included), 1L)
      to_add    <- sample(which(!included), 1L)

      proposed <- theta
      proposed[((to_add    - 1L)*dpg + 1L):(to_add    * dpg)] <-
        theta[((to_remove - 1L)*dpg + 1L):(to_remove * dpg)]
      proposed[((to_remove - 1L)*dpg + 1L):(to_remove * dpg)] <- 0

      proposed_lpl <- private$log_pl(proposed)

      # Prior and penalty cancel (same number of active entries with same values)
      logA <- proposed_lpl - lpl
      if (is.nan(logA)) logA <- -Inf

      if (log(runif(1)) < logA) {
        private$.included[to_remove] <- FALSE
        private$.included[to_add]    <- TRUE
        list(theta = proposed, lpl = proposed_lpl)
      } else {
        list(theta = theta, lpl = lpl)
      }
    },

    # ------------------------------------------------------------------
    # Move: jump — toggle a randomly chosen position
    # ------------------------------------------------------------------
    step_jump = function(theta, lpl) {
      dpg      <- private$.dim_per_group
      included <- private$.included
      npos     <- private$.npos

      jump_pos <- sample(npos, 1L)
      vec_jump <- seq_len(npos) == jump_pos
      jump_idx <- rep(vec_jump, each = dpg)

      if (included[jump_pos]) {
        # ---- Delete the selected position ----
        proposed     <- theta * as.numeric(rep(included & !vec_jump, each = dpg))
        proposed_lpl <- private$log_pl(proposed)

        logA <- (proposed_lpl - lpl) +
          sum(dnorm(proposed[proposed != 0], sd = private$.sdprior, log = TRUE)) -
          sum(dnorm(theta[theta != 0],       sd = private$.sdprior, log = TRUE)) +
          sum(dnorm(theta[jump_idx],         sd = private$.sdbirth, log = TRUE)) +
          private$.logpenalty

        if (log(runif(1)) < logA) {
          private$.included[jump_pos] <- FALSE
          list(theta = proposed, lpl = proposed_lpl)
        } else {
          list(theta = theta, lpl = lpl)
        }
      } else {
        # ---- Add the selected position ----
        proposed     <- theta + rnorm(private$.fdim, 0, private$.sdbirth) * jump_idx
        proposed_lpl <- private$log_pl(proposed)

        logA <- (proposed_lpl - lpl) +
          sum(dnorm(proposed[proposed != 0], sd = private$.sdprior, log = TRUE)) -
          sum(dnorm(theta[theta != 0],       sd = private$.sdprior, log = TRUE)) -
          sum(dnorm(proposed[jump_idx],      sd = private$.sdbirth, log = TRUE)) -
          private$.logpenalty

        if (log(runif(1)) < logA) {
          private$.included[jump_pos] <- TRUE
          list(theta = proposed, lpl = proposed_lpl)
        } else {
          list(theta = theta, lpl = lpl)
        }
      }
    }
  ),

  active = list(

    #' @field included
    #' Current logical vector of length `npos` indicating which positions in the
    #' maximal `mrfi` are active. Read-only.
    included = function() private$.included,

    #' @field n_included Number of currently active positions.
    n_included = function() sum(private$.included),

    #' @field logpenalty The log-penalty per active position (read-only).
    logpenalty = function() private$.logpenalty,

    #' @field sdbirth SD of the birth proposal distribution (read-only).
    sdbirth = function() private$.sdbirth,

    #' @field samples
    #' A [`tibble::tibble`] with columns `t`, `position`, `interaction`, and
    #' `value` containing only the **non-zero** theta entries across all
    #' collected iterations (i.e., only entries for active positions at each
    #' step). Returns `NULL` if `$run()` has not been called yet.
    samples = function() {
      if (nrow(private$.chain) == 0L) return(NULL)
      resdf <- as.data.frame(private$.chain)
      resdf$t <- private$.t_chain
      resdf <- tidyr::pivot_longer(resdf, cols = -"t")
      desc  <- mrf2d::vec_description(private$.mrfi, private$.family,
                                      private$.C)
      resdf <- cbind(resdf, desc)
      resdf <- tibble::as_tibble(resdf[, c("t", "position", "interaction",
                                           "value")])
      resdf[resdf$value != 0, ]
    },

    #' @field samples_mrfi
    #' A [`tibble::tibble`] with columns `t`, `position`, and `value`
    #' (logical) recording the active set at every collected iteration.
    #' Returns `NULL` if `$run()` has not been called yet.
    samples_mrfi = function() {
      if (nrow(private$.included_chain) == 0L) return(NULL)
      nsamples   <- nrow(private$.included_chain)
      desc       <- mrf2d::vec_description(private$.mrfi, private$.family,
                                           private$.C)
      pos_labels <- as.character(
        desc$position[seq(1L, private$.fdim, by = private$.dim_per_group)]
      )
      df         <- as.data.frame(private$.included_chain)
      colnames(df) <- pos_labels
      df$t       <- private$.t_chain
      df         <- tidyr::pivot_longer(df,
                                        cols     = tidyr::all_of(pos_labels),
                                        names_to = "position")
      tibble::as_tibble(df[, c("t", "position", "value")])
    }
  ),

  public = list(

    #' @description
    #' Create a new `MRFRJPseudoBayes` sampler.
    #'
    #' @param z A matrix with integer values in `0, ..., C` (the observed
    #'   random field).
    #' @param mrfi An `mrfi` object specifying the **maximal** interaction
    #'   structure (the algorithm can explore any subset of its positions).
    #' @param family A character string naming the parameter restriction family.
    #' @param sdprior Standard deviation of the Gaussian prior on each active
    #'   `theta` entry. Defaults to `1`.
    #' @param sdkernel Standard deviation of the `within`-move proposal.
    #'   Defaults to `0.005`.
    #' @param sdbirth Standard deviation of the `birth` and `jump`-add
    #'   proposals. Defaults to `0.05`.
    #' @param logpenalty Log-penalty per active position in the BIC-like
    #'   term. Defaults to `log(nrow(z) * ncol(z))`.
    #' @param kernel_probs Numeric vector of length 5 with unnormalised weights
    #'   for the moves `(within, swap, death, birth, jump)`. Defaults to
    #'   `c(4, 1, 1, 1, 1)`.
    #' @param init_included Initial active set. Either `"zero"` (no positions
    #'   active, default), `"full"` (all positions active), `"nearest"` (nearest neighbors)
    #'   or a logical vector of length equal to the number of positions in `mrfi`.
    #' @param init_theta Initial `theta` vector. Either:
    #'   * `"zero"` — active positions initialised with `N(0, sdkernel^2)`;
    #'   * `"pl"` — initialised at the maximum pseudolikelihood estimate
    #'     (restricted to the initially active set);
    #'   * a numeric vector of length `fdim`.
    initialize = function(z, mrfi, family,
                          sdprior       = 1,
                          sdkernel      = 0.005,
                          sdbirth       = 0.05,
                          logpenalty    = NULL,
                          kernel_probs  = c(4, 1, 1, 1, 1),
                          init_included = "zero",
                          init_theta    = "zero") {
      # Parent validates z, sets .z, .mrfi, .family, .C, .fdim
      # Use init_theta = "zero" here; we override theta below
      super$initialize(z, mrfi, family,
                       sdprior   = sdprior,
                       sdkernel  = sdkernel,
                       init_theta = "zero")

      npos <- length(mrfi)
      stopifnot(private$.fdim %% npos == 0L)
      dpg <- as.integer(private$.fdim / npos)

      # ---- Resolve included -----------------------------------------------
      if (is.character(init_included)) {
        if (init_included == "zero") {
          included <- rep(FALSE, npos)
        } else if (init_included == "full") {
          included <- rep(TRUE, npos)
        } else if(init_included == "nearest"){
          included <- sapply(as.list(mrfi), function(x) all(c(0,1) %in% x))
        } else {
          stop("'init_included' must be \"zero\", \"full\", or a logical vector.")
        }
      } else {
        stopifnot(is.logical(init_included), length(init_included) == npos)
        included <- init_included
      }

      # ---- Resolve initial theta ------------------------------------------
      if (is.character(init_theta)) {
        if (init_theta == "zero") {
          theta <- rnorm(private$.fdim, 0, sdkernel) *
            as.numeric(rep(included, each = dpg))
        } else if (init_theta == "pl") {
          theta <- mrf2d::smr_array(mrf2d::fit_pl(z, mrfi, family)$theta,
                                    family)
          theta <- theta * as.numeric(rep(included, each = dpg))
        } else {
          stop("'init_theta' must be \"zero\", \"pl\", or a numeric vector.")
        }
      } else {
        stopifnot(is.numeric(init_theta), length(init_theta) == private$.fdim)
        theta <- init_theta * as.numeric(rep(included, each = dpg))
      }
      private$.theta <- theta

      # ---- Store RJ-specific fields ---------------------------------------
      private$.included       <- included
      private$.npos           <- npos
      private$.dim_per_group  <- dpg
      private$.logpenalty     <- if (is.null(logpenalty)) log(prod(dim(z)))
                                 else as.numeric(logpenalty)
      private$.sdbirth        <- sdbirth
      private$.kernel_probs   <- kernel_probs
      private$.included_chain <- matrix(logical(0), nrow = 0L, ncol = npos)
      private$.t_chain        <- integer(0)
    },

    #' @description
    #' Run the Reversible-Jump sampler.
    #'
    #' At each iteration one move type is sampled according to `kernel_probs`
    #' and the corresponding proposal is accepted or rejected via the
    #' Metropolis-Hastings criterion. Samples are appended to previously
    #' collected samples so this method can be called repeatedly.
    #'
    #' If `warmup > 0`, the sampler first runs `warmup` iterations in which
    #' only the `within` move is proposed (i.e. `kernel_probs = c(1,0,0,0,0)`),
    #' letting the active `theta` entries settle before the interaction
    #' structure is allowed to change. These warmup iterations are stored
    #' with negative iteration indices and are excluded from inclusion
    #' probabilities and trace segmentation.
    #'
    #' @param nsamples Number of (post-warmup) iterations to run.
    #' @param warmup Number of `within`-only warmup iterations to run before
    #'   `nsamples` regular iterations. Defaults to `0`.
    #' @param verbose If `TRUE`, prints iteration progress.
    #'
    #' @return The sampler itself, invisibly (allows chaining).
    run = function(nsamples, warmup = 0, verbose = interactive()) {
      stopifnot(is.numeric(nsamples), length(nsamples) == 1L, nsamples >= 1L)
      stopifnot(is.numeric(warmup), length(warmup) == 1L, warmup >= 0L)

      total              <- warmup + nsamples
      new_theta_chain    <- matrix(0,       nrow = total, ncol = private$.fdim)
      new_included_chain <- matrix(FALSE,   nrow = total, ncol = private$.npos)

      move_list    <- c("within", "swap", "death", "birth", "jump")

      theta <- private$.theta
      lpl   <- private$log_pl(theta)

      if (verbose)
        pb <- cli::cli_progress_bar(
          total       = total,
          format      = "{cli::pb_bar} {cli::pb_current}/{cli::pb_total} | {cli::pb_rate} | ETA: {cli::pb_eta}",
          .auto_close = FALSE
        )

      for (i in seq_len(total)) {
        move <- if (i <= warmup) "within"
                else sample(move_list, 1L, prob = private$.kernel_probs)

        res <- switch(move,
          within = private$step_within(theta, lpl),
          birth  = private$step_birth (theta, lpl),
          death  = private$step_death (theta, lpl),
          swap   = private$step_swap  (theta, lpl),
          jump   = private$step_jump  (theta, lpl)
        )

        theta <- res$theta
        lpl   <- res$lpl

        new_theta_chain[i, ]    <- theta
        new_included_chain[i, ] <- private$.included

        if (verbose) cli::cli_progress_update(id = pb)
      }
      if (verbose) cli::cli_progress_done(id = pb)

      private$.theta          <- theta
      private$.chain          <- rbind(private$.chain, new_theta_chain)
      private$.included_chain <- rbind(private$.included_chain,
                                       new_included_chain)
      private$.t_chain        <- c(private$.t_chain,
                                   private$next_t(warmup, nsamples))

      invisible(self)
    },

    #' @description
    #' Compute posterior summary statistics.
    #'
    #' Returns one row per `(position, interaction)` parameter. The `prob`
    #' column is the posterior **inclusion probability** of the position
    #' (fraction of post-burnin iterations where the position was active),
    #' computed from the full inclusion chain so that positions with
    #' `prob = 0` also appear. The columns `q025`, `mean`, `q975`, and `sd`
    #' are **conditional** posterior statistics of the non-zero `theta`
    #' entries; they are `NA` for positions that were never active.
    #'
    #' @param burnin Fraction (if `< 1`) or absolute number of initial samples
    #'   to discard. Defaults to `0.25`.
    #'
    #' @return A [`tibble::tibble`] with columns `position`, `interaction`,
    #'   `prob`, `q025`, `mean`, `q975`, and `sd`.
    summary = function(burnin = 0.25) {
      if (nrow(private$.chain) == 0L)
        stop("No samples available. Call $run() first.")
      tmax   <- sum(private$.t_chain > 0L)  # warmup iterations never count
      if (burnin < 1) burnin <- floor(burnin * tmax)
      n_post <- tmax - burnin
      stopifnot(n_post > 0)

      # --- Inclusion probability from the full chain (all positions) --------
      im        <- self$samples_mrfi
      im        <- im[im$t > 0L & im$t > burnin, ]
      incl_prob <- dplyr::summarize(
        dplyr::group_by(im, .data$position),
        prob = mean(.data$value),
        .groups = "drop"
      )

      # --- Conditional theta statistics (active entries only) ---------------
      smp <- self$samples
      if (!is.null(smp)) smp <- smp[smp$t > 0L & smp$t > burnin, ]
      if (!is.null(smp) && nrow(smp) > 0L) {
        theta_stats <- dplyr::summarize(
          dplyr::group_by(smp, .data$position, .data$interaction),
          q025 = quantile(.data$value, probs = 0.025),
          mean = mean(.data$value),
          q975 = quantile(.data$value, probs = 0.975),
          sd   = sd(.data$value),
          .groups = "drop"
        )
      } else {
        theta_stats <- tibble::tibble(
          position    = character(0L), interaction = character(0L),
          q025 = numeric(0L), mean = numeric(0L),
          q975 = numeric(0L), sd   = numeric(0L)
        )
      }

      # --- Full (position, interaction) grid --------------------------------
      desc     <- mrf2d::vec_description(private$.mrfi, private$.family,
                                         private$.C)
      all_desc <- tibble::tibble(
        position    = as.character(desc$position),
        interaction = as.character(desc$interaction)
      )

      # Join: all pairs + inclusion prob + conditional theta stats
      result <- dplyr::left_join(all_desc, incl_prob, by = "position")
      dplyr::left_join(result, theta_stats,
                       by = c("position", "interaction"))
    },

    #' @description
    #' Plot the neighbourhood coloured by posterior inclusion probability.
    #'
    #' Each candidate position from the maximal `mrfi` is shown as a square
    #' at its relative offset `(rx, ry)` from the reference pixel. Squares
    #' are coloured from white (inclusion probability = 0) to red
    #' (probability = 1). Both a position and its symmetric counterpart
    #' `(-rx, -ry)` share the same colour. The reference position `(0, 0)`
    #' is drawn in black.
    #'
    #' @param burnin Fraction or absolute number of initial samples to
    #'   discard. Defaults to `0.25`.
    #'
    #' @return A `ggplot` object.
    plot_mrfi = function(burnin = 0.25) {
      if (nrow(private$.chain) == 0L)
        stop("No samples available. Call $run() first.")
      tmax   <- sum(private$.t_chain > 0L)  # warmup iterations never count
      if (burnin < 1) burnin <- floor(burnin * tmax)

      # Inclusion probabilities per position
      im        <- self$samples_mrfi
      im        <- im[im$t > 0L & im$t > burnin, ]
      incl_prob <- dplyr::summarize(
        dplyr::group_by(im, .data$position),
        prob = mean(.data$value),
        .groups = "drop"
      )

      # Parse position strings "(rx,ry)" -> integer coordinates
      coords <- t(vapply(incl_prob$position, function(s) {
        as.integer(strsplit(gsub("[() ]", "", s), ",")[[1L]])
      }, integer(2L)))

      pos_df <- tibble::tibble(
        rx   = coords[, 1L],
        ry   = coords[, 2L],
        prob = incl_prob$prob
      )
      # Add symmetric counterparts (-rx, -ry) with the same probability
      sym_df  <- tibble::tibble(rx = -pos_df$rx, ry = -pos_df$ry,
                                prob = pos_df$prob)
      tile_df <- dplyr::distinct(rbind(pos_df, sym_df))

      center_df <- tibble::tibble(rx = 0L, ry = 0L)

      ggplot2::ggplot() +
        ggplot2::geom_tile(
          data    = tile_df,
          mapping = ggplot2::aes(x = .data$rx, y = .data$ry,
                                 fill = .data$prob),
          color   = "grey50", width = 0.9, height = 0.9
        ) +
        ggplot2::scale_fill_gradient(
          name   = "P(included)",
          low    = "white", high = "red",
          limits = c(0, 1)
        ) +
        ggplot2::geom_tile(
          data    = center_df,
          mapping = ggplot2::aes(x = .data$ry, y = .data$rx),
          fill    = "black", color = "grey50", width = 0.9, height = 0.9
        ) +
        ggplot2::coord_equal() +
        ggplot2::labs(x = "ry", y = "rx") +
        ggplot2::theme_bw()
    },

    #' @description
    #' Plot the MCMC traces of the active `theta` entries.
    #'
    #' Line segments are broken wherever a position leaves and re-enters the
    #' active set, so runs of contiguous activity appear as separate segments.
    #' Iterations are thinned before plotting.
    #'
    #' @param burnin Fraction or absolute number of initial samples to exclude
    #'   from credible-band computations. Defaults to `0`.
    #' @param thin Keep every `thin`-th iteration for the plot. Defaults to `50`.
    #'
    #' @return A `ggplot` object.
    plot = function(burnin = 0, thin = 50L) {
      df <- self$samples
      if (is.null(df) || nrow(df) == 0L)
        stop("No samples available. Call $run() first.")

      df <- dplyr::group_by(df, .data$position, .data$interaction)
      df <- dplyr::mutate(df,
        dif = c(0L, diff(.data$t) - 1L),
        grp = cumsum(.data$dif)
      )
      df <- dplyr::ungroup(df)
      df <- dplyr::filter(df, .data$t %% thin == 0L)

      ggplot2::ggplot(df, ggplot2::aes(
        x = .data$t, y = .data$value, color = .data$position
      )) +
        ggplot2::geom_line(ggplot2::aes(
          group = interaction(.data$grp, .data$position)
        )) +
        ggplot2::geom_hline(yintercept = 0) +
        ggplot2::theme_bw() +
        ggplot2::facet_wrap(~.data$interaction)
    },

    #' @description Print a compact description of the sampler.
    print = function(...) {
      cat(glue::glue(
        "MRFRJPseudoBayes sampler\n",
        "  Max. interaction : {mrf2d::mrfi_to_string(private$.mrfi)}\n",
        "  Active positions : {sum(private$.included)}/{private$.npos}\n",
        "  Family           : {private$.family}\n",
        "  C                : {private$.C}\n",
        "  Samples          : {nrow(private$.chain)}\n",
        "  sdprior          : {private$.sdprior}\n",
        "  sdkernel         : {private$.sdkernel}\n",
        "  sdbirth          : {private$.sdbirth}\n",
        "  logpenalty       : {round(private$.logpenalty, 3)}\n"
      ))
      invisible(self)
    }
  )
)
