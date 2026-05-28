#' @title Metropolis-within-Gibbs sampler for Hidden MRF with Gaussian emissions
#'
#' @description
#' An R6 class extending [MRFPseudoBayes] for Bayesian inference of Hidden
#' Markov Random Fields (HMRF) with Gaussian emission distributions.
#'
#' @details
#' The model assumes a latent categorical field `z` governed by a Markov
#' Random Field with parameter array `theta`, and a continuous observed field
#' `y` whose entries follow class-specific Gaussian distributions:
#' \deqn{y_{ij} \mid z_{ij} = k \sim N(\mu_k, \sigma^2_k)}
#'
#' Conjugate priors are placed on the emission parameters:
#' \deqn{\mu_k \sim N(0, \sigma^2_{\mu}), \quad 1/\sigma^2_k \sim \text{Gamma}(\alpha_1, \alpha_2)}
#'
#' At each MCMC iteration the following updates are performed in order:
#' 1. **Emission parameters**: `mu_k` and `sigma2_k` are updated from their
#'    conjugate full conditionals given the current `z` and `y`.
#' 2. **Latent field**: `z` is updated via a single-cycle Gibbs sampler whose
#'    conditional weights combine emission log-likelihoods and the MRF prior.
#' 3. **MRF parameter**: `theta` is updated via the same Metropolis-Hastings
#'    pseudo-posterior step used in [MRFPseudoBayes] with the newly sampled `z`.
#'
#' The chain state is preserved between calls to `$run()`.
#'
#' @importFrom mrf2d expand_array smr_stat smr_array fit_pl pl_mrf2d fit_ghm
#'   vec_description mrfi_to_string
#' @importFrom tidyr pivot_longer
#' @importFrom tibble as_tibble
#' @importFrom ggplot2 ggplot aes geom_line geom_tile scale_fill_gradient
#'   theme_bw facet_wrap
#' @importFrom dplyr group_by summarize bind_rows
#' @importFrom stats dnorm rnorm runif quantile sd rgamma
#' @importFrom glue glue
#'
#' @export
HMRFPseudoBayes <- R6::R6Class(
  classname = "HMRFPseudoBayes",
  inherit = MRFPseudoBayes,

  private = list(
    .y            = NULL,  # observed continuous field (matrix)
    .C_val        = NULL,  # copy of C for use inside child private methods
    .mu           = NULL,  # current emission means, length C+1
    .sigma2       = NULL,  # current emission variances, length C+1
    .sdpriormu    = NULL,  # sd of Gaussian prior on mu
    .alpha1       = NULL,  # shape of Gamma prior on precision
    .alpha2       = NULL,  # rate  of Gamma prior on precision
    .mu_chain     = NULL,  # nsamples x (C+1) matrix
    .sigma2_chain = NULL,  # nsamples x (C+1) matrix
    .z_counts     = NULL,  # array [nrow, ncol, C+1] accumulator

    # ------------------------------------------------------------------
    # Conjugate update of mu and sigma2 given current z and y
    # ------------------------------------------------------------------
    step_pars = function() {
      C_val     <- private$.C
      y_vec     <- as.vector(private$.y)
      z_vec     <- as.vector(private$.z)
      sigma2_mu <- private$.sdpriormu^2
      sigma2_c  <- private$.sigma2

      # Sufficient statistics per label k in 0:C
      stts <- vapply(0:C_val, function(k) {
        x <- y_vec[z_vec == k]
        c(sum(x), length(x), sum(x^2))
      }, numeric(3))
      sums <- stts[1L, ]
      ns   <- stts[2L, ]
      sqs  <- stts[3L, ]

      # mu | z, sigma2, y  (Gaussian conjugate)
      mu_new <- rnorm(
        C_val + 1L,
        mean = sums / (ns + sigma2_c / sigma2_mu),
        sd   = (sqrt(1 / sigma2_mu + ns / sigma2_c))^(-1)
      )
      private$.mu <- mu_new

      # sigma2 | z, mu, y  (Inverse-Gamma conjugate)
      sqdifs  <- sqs - 2 * sums * mu_new + ns * mu_new^2
      invsigs <- rgamma(
        C_val + 1L,
        shape = private$.alpha1 + ns / 2,
        rate  = private$.alpha2 + 0.5 * sqdifs
      )
      private$.sigma2 <- 1 / invsigs
    },

    # ------------------------------------------------------------------
    # Single-cycle Gibbs update of the latent field z
    # ------------------------------------------------------------------
    step_z = function() {
      C_val  <- private$.C
      mu     <- private$.mu
      sigma2 <- private$.sigma2
      y      <- private$.y

      # Log emission weights: numeric array [nrow, ncol, C+1]
      w_list <- lapply(seq_len(C_val + 1L), function(k) {
        dnorm(y, mean = mu[k], sd = sqrt(sigma2[k]), log = TRUE)
      })
      w <- array(do.call(c, w_list), dim = c(dim(y), C_val + 1L))

      private$.z <- inner_gibbs_conditional(
        private$.z,
        cond_weights = w,
        R            = private$.mrfi@Rmat,
        ncycles      = 1L,
        theta        = mrf2d::expand_array(
          private$.theta, private$.family, private$.mrfi, C_val
        )
      )
    }
  ),

  active = list(

    #' @field y The observed continuous field (read-only).
    y = function() private$.y,

    #' @field mu Current emission means (read-only numeric vector).
    mu = function() private$.mu,

    #' @field sigma2 Current emission variances (read-only numeric vector).
    sigma2 = function() private$.sigma2,

    #' @field sdpriormu SD of the Gaussian prior on `mu` (read-only).
    sdpriormu = function() private$.sdpriormu,

    #' @field alpha1 Shape of the Gamma prior on precision (read-only).
    alpha1 = function() private$.alpha1,

    #' @field alpha2 Rate of the Gamma prior on precision (read-only).
    alpha2 = function() private$.alpha2,

    #' @field samples_pars
    #' A [`tibble::tibble`] with columns `t`, `k`, `par`, and `value`
    #' containing all collected samples of emission parameters (`mu` and
    #' `sigma2`). Returns `NULL` if `$run()` has not been called yet.
    samples_pars = function() {
      if (nrow(private$.mu_chain) == 0L) return(NULL)
      C_val    <- private$.C
      nsamples <- nrow(private$.mu_chain)
      k_labels <- as.character(0:C_val)

      dfmu       <- as.data.frame(private$.mu_chain)
      colnames(dfmu) <- k_labels
      dfmu$t   <- seq_len(nsamples)
      dfmu$par <- "mu"

      dfsig       <- as.data.frame(private$.sigma2_chain)
      colnames(dfsig) <- k_labels
      dfsig$t   <- seq_len(nsamples)
      dfsig$par <- "sigma2"

      df  <- rbind(dfmu, dfsig)
      df  <- tidyr::pivot_longer(df, cols = tidyr::all_of(k_labels), names_to = "k")
      df$k <- as.integer(df$k)
      tibble::as_tibble(df[, c("t", "k", "par", "value")])
    },

    #' @field zprobs
    #' Posterior label probability array with dimensions
    #' `c(nrow(y), ncol(y), C+1)`. The `[i, j, k+1]` entry is the empirical
    #' posterior probability that `z[i,j] == k`. Returns `NULL` before
    #' `$run()`.
    zprobs = function() {
      n <- nrow(private$.chain)
      if (n == 0L) return(NULL)
      private$.z_counts / n
    }
  ),

  public = list(

    #' @description
    #' Create a new `HMRFPseudoBayes` sampler.
    #'
    #' @param y A numeric matrix of observed continuous values.
    #' @param mrfi An `mrfi` object specifying the interaction structure for
    #'   the latent MRF.
    #' @param family A character string naming the parameter restriction
    #'   family (e.g. `"onepar"`, `"absdif"`, `"free"`).
    #' @param C Non-negative integer; the field takes values `0, ..., C`.
    #' @param sdprior Standard deviation of the Gaussian prior on each
    #'   element of `theta`. Defaults to `10`.
    #' @param sdkernel Standard deviation of the isotropic Gaussian
    #'   random-walk proposal for `theta`. Defaults to `0.005`.
    #' @param sdpriormu Standard deviation of the Gaussian prior on each
    #'   `mu_k`. Defaults to `10`.
    #' @param alpha1 Shape parameter of the Gamma prior on each precision
    #'   `1/sigma2_k`. Defaults to `0.001` (weakly informative).
    #' @param alpha2 Rate parameter of the Gamma prior on each precision.
    #'   Defaults to `0.001`.
    #' @param z0 Integer matrix; initial latent field. If `NULL`, initial
    #'   values are computed via the EM algorithm ([mrf2d::fit_ghm]).
    #' @param mu0 Numeric vector of length `C+1`; initial emission means.
    #'   Computed via EM when `NULL` and `z0` is also `NULL`.
    #' @param sigma2_0 Numeric vector of length `C+1`; initial emission
    #'   variances. Computed via EM when `NULL` and `z0` is also `NULL`.
    #' @param init_theta Initial value for `theta`. Same options as in
    #'   [MRFPseudoBayes]: `"zero"` (default), `"pl"`, or a numeric vector.
    #' @param verbose If `TRUE`, prints progress during EM initialisation.
    initialize = function(y, mrfi, family, C,
                          sdprior   = 10,
                          sdkernel  = 0.005,
                          sdpriormu = 10,
                          alpha1    = 0.001,
                          alpha2    = 0.001,
                          z0        = NULL,
                          mu0       = NULL,
                          sigma2_0  = NULL,
                          init_theta = "zero",
                          verbose    = interactive()) {
      stopifnot(is.matrix(y), is.numeric(y))
      C <- as.integer(C)
      stopifnot(C >= 1L)

      # ---- Initial values via EM if not provided --------------------------
      if (is.null(z0) || is.null(mu0) || is.null(sigma2_0)) {
        if (verbose) message("Computing initial values via EM algorithm...")
        theta_em <- array(0, dim = c(C + 1L, C + 1L, length(mrfi)))
        EM <- mrf2d::fit_ghm(
          Y          = y,
          mrfi       = mrfi,
          theta      = theta_em,
          equal_vars = TRUE,
          verbose    = FALSE
        )
        if (is.null(z0))      z0       <- EM$Z_pred
        if (is.null(mu0))     mu0      <- EM$par$mu
        if (is.null(sigma2_0)) sigma2_0 <- EM$par$sigma^2
        if (verbose) message("Done!")
      }

      stopifnot(is.matrix(z0), length(mu0) == C + 1L,
                length(sigma2_0) == C + 1L, all(sigma2_0 > 0))

      # ---- Initialise theta/z chain via parent ----------------------------
      super$initialize(z0, mrfi, family,
                       sdprior    = sdprior,
                       sdkernel   = sdkernel,
                       init_theta = init_theta)

      private$.y            <- y
      private$.mu           <- mu0
      private$.sigma2       <- sigma2_0
      private$.sdpriormu    <- sdpriormu
      private$.alpha1       <- alpha1
      private$.alpha2       <- alpha2
      private$.mu_chain     <- matrix(numeric(0), nrow = 0L, ncol = C + 1L)
      private$.sigma2_chain <- matrix(numeric(0), nrow = 0L, ncol = C + 1L)
      private$.z_counts     <- array(0, dim = c(dim(y), C + 1L))
    },

    #' @description
    #' Run the Metropolis-within-Gibbs sampler.
    #'
    #' Each iteration performs in order: (1) conjugate update of emission
    #' parameters, (2) Gibbs update of the latent field, (3) MH update of
    #' `theta` using the newly sampled field. Samples are appended to any
    #' previously collected samples.
    #'
    #' @param nsamples Number of iterations to run.
    #' @param verbose If `TRUE`, prints iteration progress to the console.
    #'
    #' @return The sampler itself, invisibly (allows chaining).
    run = function(nsamples, verbose = interactive()) {
      stopifnot(is.numeric(nsamples), length(nsamples) == 1L, nsamples >= 1L)

      C_val            <- private$.C
      new_theta_chain  <- matrix(0, nrow = nsamples, ncol = private$.fdim)
      new_mu_chain     <- matrix(0, nrow = nsamples, ncol = C_val + 1L)
      new_sigma2_chain <- matrix(0, nrow = nsamples, ncol = C_val + 1L)

      theta         <- private$.theta
      log_target_val <- private$log_target(theta)

      for (i in seq_len(nsamples)) {
        # 1. Update emission parameters (mu, sigma2)
        private$step_pars()

        # 2. Update latent field z
        private$step_z()
        private$.z_counts <- private$.z_counts +
          indicator_array(private$.z, C_val)

        # 3. Update theta via MH (log_target re-evaluated with new z)
        log_target_val <- private$log_target(theta)
        result         <- private$mh_step(theta, log_target_val)
        theta          <- result$theta
        log_target_val <- result$log_target

        new_mu_chain[i, ]     <- private$.mu
        new_sigma2_chain[i, ] <- private$.sigma2
        new_theta_chain[i, ]  <- theta

        if (verbose) cat("\r", nrow(private$.chain) + i)
      }
      if (verbose) cat("\n")

      private$.theta        <- theta
      private$.chain        <- rbind(private$.chain, new_theta_chain)
      private$.mu_chain     <- rbind(private$.mu_chain,     new_mu_chain)
      private$.sigma2_chain <- rbind(private$.sigma2_chain, new_sigma2_chain)

      invisible(self)
    },

    #' @description
    #' Compute posterior summary statistics for the emission parameters.
    #'
    #' @param burnin Fraction (if `< 1`) or absolute number of initial
    #'   samples to discard. Defaults to `0.25` (25% burn-in).
    #'
    #' @return A [`tibble::tibble`] with columns `k`, `par`, `q025`, `mean`,
    #'   `q975`, and `sd`.
    summary_pars = function(burnin = 0.25) {
      df <- self$samples_pars
      if (is.null(df)) stop("No samples available. Call $run() first.")
      tmax <- max(df$t)
      if (burnin < 1) burnin <- burnin * tmax
      df <- df[df$t > burnin, ]
      dplyr::summarize(
        dplyr::group_by(df, .data$k, .data$par),
        q025 = quantile(.data$value, probs = 0.025),
        mean = mean(.data$value),
        q975 = quantile(.data$value, probs = 0.975),
        sd   = sd(.data$value),
        .groups = "drop"
      )
    },

    #' @description
    #' Plot MCMC traces or posterior label probability maps.
    #'
    #' @param what One of:
    #'   * `"pars"` — trace plots of `mu` and `sigma2` (default);
    #'   * `"theta"` — trace plot of `theta` (delegates to [MRFPseudoBayes]);
    #'   * `"zprobs"` — heatmaps of posterior label probabilities.
    #' @param burnin Burn-in passed to `$summary()` when `what = "theta"`.
    #'
    #' @return A `ggplot` object.
    plot = function(what = "pars", burnin = 0) {
      if (what == "theta") return(super$plot(burnin = burnin))

      if (what == "pars") {
        df <- self$samples_pars
        if (is.null(df)) stop("No samples available. Call $run() first.")
        df$k <- as.factor(df$k)
        return(
          ggplot2::ggplot(df, ggplot2::aes(
            x = .data$t, y = .data$value, color = .data$k
          )) +
            ggplot2::geom_line() +
            ggplot2::facet_wrap(~par, scales = "free") +
            ggplot2::theme_bw()
        )
      }

      if (what == "zprobs") {
        zp <- self$zprobs
        if (is.null(zp)) stop("No samples available. Call $run() first.")
        C_val  <- private$.C
        slices <- lapply(seq_len(C_val + 1L), function(k) {
          sl <- zp[, , k]
          data.frame(
            x     = as.vector(row(sl)),
            y     = as.vector(col(sl)),
            value = as.vector(sl),
            k     = k - 1L
          )
        })
        df <- dplyr::bind_rows(slices)
        return(
          ggplot2::ggplot(df, ggplot2::aes(
            x = .data$x, y = .data$y, fill = .data$value
          )) +
            ggplot2::geom_tile() +
            ggplot2::scale_fill_gradient(low = "white", high = "black") +
            ggplot2::facet_wrap(~k) +
            ggplot2::theme_bw()
        )
      }

      stop("'what' must be one of \"pars\", \"theta\", or \"zprobs\".")
    },

    #' @description Print a compact description of the sampler.
    print = function(...) {
      cat(glue::glue(
        "HMRFPseudoBayes sampler\n",
        "  Interaction  : {mrf2d::mrfi_to_string(private$.mrfi)}\n",
        "  Family       : {private$.family}\n",
        "  C            : {private$.C}\n",
        "  Field size   : {nrow(private$.y)} x {ncol(private$.y)}\n",
        "  Samples      : {nrow(private$.chain)}\n",
        "  sdprior      : {private$.sdprior}\n",
        "  sdkernel     : {private$.sdkernel}\n",
        "  sdpriormu    : {private$.sdpriormu}\n"
      ))
      invisible(self)
    }
  )
)
