#' @param An observed continuous-valued `numeric matrix` of values.
#'
#' @return A `list` with ...
#' @export
hmrfbayes <- function(y, llapprox, nsamples = 1000, 
                      z0 = NULL, mu0 = NULL, sigma0 = NULL,
                      theta0 = NULL, sdpriortheta = 1,
                      sdkerneltheta = 0.005, sdpriormu = 10,
                      alpha1 = 0.001, alpha2 = 0.001,
                      verbose = interactive()){
  C <- llapprox@C
  mrfi <- llapprox@mrfi
  family <- llapprox@family

  # Check and set initial field config if not defined
  if(is.null(z0) | is.null(mu0) | is.null(sigma0) | is.null(theta0)){
    if(verbose) cat("Computing initial values via EM algorithm", "\n")
    EM <- mrf2d::fit_ghm(Y = y, 
                         mrfi = mrfi(1) - c(1,0),
                         theta = array(0, dim = c(C+1, C+1, 1)),
                         equal_vars = TRUE,
                         verbose = FALSE)
    if(is.null(z0)) z0 <- EM$Z_pred
    if(is.null(mu0)) mu0 <- EM$par$mu 
    if(is.null(sigma0)) sigma0 <- EM$par$sigma^2
    if(is.null(theta0)) 
      theta0 <- mrf2d::smr_array(mrf2d::fit_pl(EM$Z_pred, mrfi, family)$theta,
                                 family)
    if(verbose) cat("Done!", "\n")
  }

  if(is.array(theta0)) theta0 <- mrf2d::smr_array(theta0, family)
  dimtheta <- length(theta0)
  
  # Initialize structures to store the results
  # Testing with a fixed theta first, this will break with
  # interactions other than nearest-neighbor.
  mu_chain <- sigma2_chain <- matrix(0, nrow = nsamples, ncol = C+1)
  theta_chain <- matrix(0, nrow = nsamples, ncol = dimtheta)
  z_current <- z0
  sigma2_mu <- sdpriormu^2

  mu_current <- mu0
  sigma2_current <- sigma0
  theta_current <- theta0
  z_counts <- array(0, c(dim(y), C+1))

  # Generate Markov Chain
  if(verbose) cat("Running Markov Chain iterations...", "\n")
  for(t in seq_len(nsamples)){
      # Update mus
      stts <- tapply(as.vector(y), INDEX = as.vector(z_current), 
          FUN = function(x){c(sum(x), length(x), sum(x^2))}) 
      sums <- sapply(stts, "[", 1)
      ns <- sapply(stts, "[", 2)
      mu_current <- rnorm(C+1,
          mean = sums/(ns + sigma2_current/sigma2_mu),
          sd = sqrt(1/sigma2_mu + ns/sigma2_current)^-1)
      mu_chain[t,] <- mu_current

      # Update Sigmas
      sqs <- sapply(stts, "[", 3)
      sqdifs <- sqs -2*sums*mu_current + ns*mu_current^2
      invsigs <- rgamma(C+1,
          shape = alpha1 + ns/2,
          rate = alpha2 + 1/2*sqdifs)
      sigma2_current <- 1/invsigs
      sigma2_chain[t,] <- sigma2_current

      # Update latent field
      w2 <- lapply(seq_len(C+1), FUN = function(x){
          dnorm(y, 
                mean = mu_current[x],
                sd = sqrt(sigma2_current[x]),
                log = TRUE)
      })
      w <- array(do.call(c, w2), dim = c(dim(y), length(mu_current)))
      z_current <- inner_gibbs_conditional(z_current,
          cond_weights = w, R = mrfi@Rmat, ncycles = 1,
          theta = mrf2d::expand_array(theta_current, family, mrfi, C))
      z_counts <- z_counts + indicator_array(z_current, C)
      # Update theta with the newly sampled field
      theta_prop <- theta_current + rnorm(dimtheta, sd = sdkerneltheta)
      if(!llapprox@pass_entire) {
        zpass <- mrf2d::smr_stat(z_current, mrfi, family)
      } else {
        zpass <- z_current
      }
      logratio <- llapprox@lafn(zpass, theta_prop) -
        llapprox@lafn(zpass, theta_current) +
        sum(dnorm(theta_prop, sd = sdpriortheta, log = TRUE)) - 
        sum(dnorm(theta_current, sd = sdpriortheta, log = TRUE))
      if(log(runif(1)) < logratio) theta_current <- theta_prop
      theta_chain[t,] <- theta_current

      if(verbose) cat("\r", t)
  }
  if(verbose) cat("\n", "Done!", "\n")

  # Format data.frame with pars
  dfmus <- as.data.frame(mu_chain)
  colnames(dfmus) <- 0:C
  dfmus <- cbind(t = 1:nsamples, dfmus, par = "mu")
  
  dfsigma <- as.data.frame(sigma2_chain)
  colnames(dfsigma) <- 0:C
  dfsigma <- cbind(t = 1:nsamples, dfsigma, par = "sigma2")

  dfpars <- rbind(dfmus, dfsigma)
  
  # Format data.frame with theta
  dftheta <- cbind(t = 1:nsamples, as.data.frame(theta_chain), 
                   mrf2d::vec_description(mrfi, family, C))
  if(family != "onepar"){
    dftheta <- tidyr::pivot_longer(dftheta, cols = -"t")
  } else {
    colnames(dftheta)[2] <- "value"
  }

  out <- list(dfpars = tibble::as_tibble(dfpars),
              dftheta = tibble::as_tibble(dftheta),
              ll = llapprox,
              last_z = z_current,
              zprobs = z_counts/nsamples)

  class(out) <- "hmrfbayes_out"
  return(out)
}
