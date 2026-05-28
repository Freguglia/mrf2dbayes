#' @title Metropolis-Hastings sampler for MRF with fixed interaction structure
#'
#' @description
#' An R6 class implementing a Metropolis-Hastings algorithm for Bayesian
#' inference of the parameter array (theta) of a Markov Random Field with a
#' fixed and known interaction structure (`mrfi`). The posterior is approximated
#' via pseudolikelihood combined with independent Gaussian priors on each
#' parameter.
#'
#' @details
#' The log pseudo-posterior used as the target is:
#' \deqn{\log \tilde{\pi}(\theta | z) = \log PL(z | \theta) + \sum_k \log \phi(\theta_k; 0, \sigma_{\text{prior}})}
#' where \eqn{PL} is the pseudolikelihood of the MRF and \eqn{\phi} is the
#' normal density.
#'
#' Proposals are generated via an isotropic normal random walk:
#' \deqn{\theta^* = \theta + \varepsilon, \quad \varepsilon \sim N(0, \sigma_{\text{kernel}}^2 I)}
#'
#' The chain state is preserved between calls to `$run()`, so the sampler can
#' be continued by calling `$run()` again.
#'
#' @importFrom mrf2d expand_array smr_stat smr_array fit_pl pl_mrf2d vec_description mrfi_to_string
#' @importFrom tidyr pivot_longer
#' @importFrom tibble as_tibble
#' @importFrom ggplot2 ggplot aes geom_line geom_rect theme_bw facet_wrap
#' @importFrom dplyr group_by summarize
#' @importFrom stats dnorm rnorm runif quantile sd
#' @importFrom glue glue
#' @importFrom cli cli_progress_bar cli_progress_update cli_progress_done
#'
#' @export
MRFPseudoBayes <- R6::R6Class(
  classname = "MRFPseudoBayes",

  private = list(
    .z        = NULL,  # observed field (matrix)
    .mrfi     = NULL,  # interaction structure (mrfi object)
    .family   = NULL,  # parameter restriction family (character)
    .C        = NULL,  # max label value (integer: field takes values 0..C)
    .fdim     = NULL,  # length of the vectorised theta
    .sdprior  = NULL,  # sd of the Gaussian prior
    .sdkernel = NULL,  # sd of the random-walk proposal kernel
    .theta    = NULL,  # current chain state (numeric vector)
    .chain    = NULL,  # collected samples (nsamples x fdim matrix)

    log_pl = function(theta_vec) {
      theta_arr <- mrf2d::expand_array(theta_vec, private$.family,
                                       private$.mrfi, private$.C)
      mrf2d::pl_mrf2d(private$.z, private$.mrfi, theta_arr, log_scale = TRUE)
    },

    log_prior = function(theta_vec) {
      sum(dnorm(theta_vec, sd = private$.sdprior, log = TRUE))
    },

    log_target = function(theta_vec) {
      private$log_pl(theta_vec) + private$log_prior(theta_vec)
    },

    # Single Metropolis-Hastings step on theta.
    # Returns a list(theta, log_target) with the (possibly updated) state.
    mh_step = function(theta, log_target_val) {
      proposed       <- theta + rnorm(private$.fdim, mean = 0,
                                      sd = private$.sdkernel)
      log_target_new <- private$log_target(proposed)
      if (log(runif(1)) < log_target_new - log_target_val) {
        list(theta = proposed, log_target = log_target_new)
      } else {
        list(theta = theta, log_target = log_target_val)
      }
    }
  ),

  active = list(

    #' @field samples
    #' A [`tibble::tibble`] with columns `t`, `position`, `interaction`, and
    #' `value` containing all collected posterior samples. Returns `NULL` if
    #' `$run()` has not been called yet.
    samples = function() {
      if (nrow(private$.chain) == 0L) return(NULL)
      resdf <- as.data.frame(private$.chain)
      resdf$t <- seq_len(nrow(resdf))
      resdf <- tidyr::pivot_longer(resdf, cols = -"t")
      desc  <- mrf2d::vec_description(private$.mrfi, private$.family,
                                      private$.C)
      resdf <- cbind(resdf, desc)
      tibble::as_tibble(resdf[, c("t", "position", "interaction", "value")])
    },

    #' @field theta Current state of the Markov chain (numeric vector).
    theta = function() private$.theta,

    #' @field n_samples Number of posterior samples collected so far.
    n_samples = function() nrow(private$.chain),

    #' @field z The observed random field. Can be replaced with a field of the
    #'   same dimensions.
    z = function(value) {
      if (missing(value)) return(private$.z)
      if (!is.matrix(value)) stop("'z' must be a matrix.")
      if (!identical(dim(value), dim(private$.z)))
        stop("Dimensions of new 'z' must match the original (",
             paste(dim(private$.z), collapse = "x"), ").")
      private$.z <- value
    },

    #' @field mrfi The interaction structure.
    mrfi = function() private$.mrfi,

    #' @field family The parameter restriction family.
    family = function() private$.family,

    #' @field sdprior SD of the Gaussian prior (read-only).
    sdprior = function() private$.sdprior,

    #' @field sdkernel SD of the random-walk proposal kernel (read-only).
    sdkernel = function() private$.sdkernel
  ),

  public = list(

    #' @description
    #' Create a new `MRFPseudoBayes` sampler.
    #'
    #' @param z A matrix with integer values in `0, ..., C` representing the
    #'   observed random field.
    #' @param mrfi An `mrfi` object specifying the (fixed) interaction structure.
    #' @param family A character string naming the parameter restriction family
    #'   (e.g. `"onepar"`, `"oneeach"`, `"absdif"`, `"dif"`, `"free"`).
    #' @param sdprior Standard deviation of the independent Gaussian prior placed
    #'   on each parameter. Defaults to `10` (weakly informative).
    #' @param sdkernel Standard deviation of the isotropic Gaussian random-walk
    #'   proposal. Defaults to `0.05`.
    #' @param init_theta Initial value of theta. Either:
    #'   * `"zero"` — start at a vector of zeros (default);
    #'   * `"pl"` — start at the maximum pseudolikelihood estimate;
    #'   * a numeric vector of length equal to the vectorised theta.
    initialize = function(z, mrfi, family,
                          sdprior   = 10,
                          sdkernel  = 0.05,
                          init_theta = "zero") {
      stopifnot(is.matrix(z))
      labels <- sort(unique(as.vector(z)))
      if (!identical(labels, seq(0L, max(labels), 1L))) {
        stop("'z' must contain integer values in 0, 1, ..., C for some C.")
      }

      private$.z      <- z
      private$.mrfi   <- mrfi
      private$.family <- family
      private$.C      <- max(labels)

      T_z <- mrf2d::smr_stat(z, mrfi, family)
      private$.fdim <- length(T_z)

      private$.sdprior  <- sdprior
      private$.sdkernel <- sdkernel

      if (is.character(init_theta)) {
        if (init_theta == "zero") {
          private$.theta <- T_z * 0
        } else if (init_theta == "pl") {
          private$.theta <- mrf2d::smr_array(
            mrf2d::fit_pl(z, mrfi, family)$theta, family
          )
        } else {
          stop("'init_theta' must be \"zero\", \"pl\", or a numeric vector.")
        }
      } else {
        stopifnot(is.numeric(init_theta), length(init_theta) == private$.fdim)
        private$.theta <- init_theta
      }

      private$.chain <- matrix(numeric(0), nrow = 0L, ncol = private$.fdim)
    },

    #' @description
    #' Run the Metropolis-Hastings sampler.
    #'
    #' Samples are appended to any previously collected samples, so this method
    #' can be called multiple times to extend the chain.
    #'
    #' @param nsamples Number of iterations to run.
    #' @param verbose If `TRUE`, prints iteration progress to the console.
    #'
    #' @return The sampler itself, invisibly (allows chaining: `$run(1000)$plot()`).
    run = function(nsamples, verbose = interactive()) {
      stopifnot(is.numeric(nsamples), length(nsamples) == 1L, nsamples >= 1L)

      new_chain  <- matrix(0, nrow = nsamples, ncol = private$.fdim)
      theta      <- private$.theta
      log_target <- private$log_target(theta)

      if (verbose)
        pb <- cli::cli_progress_bar(
          total       = nsamples,
          format      = "{cli::pb_bar} {cli::pb_current}/{cli::pb_total} | {cli::pb_rate} | ETA: {cli::pb_eta}",
          .auto_close = FALSE
        )

      for (i in seq_len(nsamples)) {
        result     <- private$mh_step(theta, log_target)
        theta      <- result$theta
        log_target <- result$log_target

        new_chain[i, ] <- theta

        if (verbose) cli::cli_progress_update(id = pb)
      }
      if (verbose) cli::cli_progress_done(id = pb)

      private$.theta <- theta
      private$.chain <- rbind(private$.chain, new_chain)

      invisible(self)
    },

    #' @description
    #' Compute posterior summary statistics.
    #'
    #' @param burnin Fraction (if `< 1`) or absolute number of initial samples
    #'   to discard. Defaults to `0.25` (25% burn-in).
    #'
    #' @return A [`tibble::tibble`] with columns `position`, `interaction`,
    #'   `q025`, `mean`, `q975`, and `sd`.
    summary = function(burnin = 0.25) {
      df <- self$samples
      if (is.null(df)) stop("No samples available. Call $run() first.")
      tmax <- max(df$t)
      if (burnin < 1) burnin <- burnin * tmax
      df <- df[df$t > burnin, ]
      dplyr::summarize(
        dplyr::group_by(df, .data$position, .data$interaction),
        q025 = quantile(.data$value, probs = 0.025),
        mean = mean(.data$value),
        q975 = quantile(.data$value, probs = 0.975),
        sd   = sd(.data$value),
        .groups = "drop"
      )
    },

    #' @description
    #' Plot the MCMC trace with posterior credible bands.
    #'
    #' @param burnin Fraction or absolute number of initial samples to exclude
    #'   from the credible bands (not from the trace). Defaults to `0`.
    #'
    #' @return A `ggplot` object.
    plot = function(burnin = 0) {
      df   <- self$samples
      if (is.null(df)) stop("No samples available. Call $run() first.")
      stts <- self$summary(burnin = burnin)
      tmax <- max(df$t)

      ggplot2::ggplot(df) +
        ggplot2::geom_line(
          ggplot2::aes(x = .data$t, y = .data$value, color = .data$position)
        ) +
        ggplot2::geom_rect(
          data = stts,
          ggplot2::aes(
            xmin  = 0, xmax  = tmax,
            ymin  = .data$q025, ymax  = .data$q975,
            fill  = .data$position
          ),
          alpha = 0.15
        ) +
        ggplot2::theme_bw() +
        ggplot2::facet_wrap(~.data$interaction)
    },

    #' @description Print a compact description of the sampler.
    print = function(...) {
      cat(glue::glue(
        "MRFPseudoBayes sampler\n",
        "  Interaction : {mrf2d::mrfi_to_string(private$.mrfi)}\n",
        "  Family      : {private$.family}\n",
        "  C           : {private$.C}\n",
        "  Samples     : {nrow(private$.chain)}\n",
        "  sdprior     : {private$.sdprior}\n",
        "  sdkernel    : {private$.sdkernel}\n"
      ))
      invisible(self)
    }
  )
)
