#' @param An observed continuous-valued `numeric matrix` of values.
#'
#' @return A `list` with ...
#' @export
hmrfbayes <- function(y, llapprox, nsamples = 1000, 
                      z0 = NULL, mu0 = NULL, sigma0 = NULL,
                      sdpriormu = 10,
                      alpha1 = 0.001, alpha2 = 0.001,
                      verbose = interactive()){
  C <- llapprox@C
  mrfi <- llapprox@mrfi

  # Check and set initial field config if not defined
  if(is.null(z0) | is.null(mu0) | is.null(sigma0)){
    if(verbose) cat("Computing initial values via EM algorithm", "\n")
    EM <- fit_ghm(Y = y, 
                  mrfi = mrfi(1) - c(1,0),
                  theta = array(0, dim = c(C+1, C+1, 1)),
                  verbose = FALSE)
    if(is.null(z0)) z0 <- EM$Z_pred
    if(is.null(mu0)) mu0 <- EM$par$mu 
    if(is.null(sigma0)) sigma0 <- EM$par$sigma^2
    if(verbose) cat("Done!", "\n")
  }
  
  # Initialize structures to store the results
  # Testing with a fixed theta first, this will break with
  # interactions other than nearest-neighbor.
  theta <- mrf2d::theta_potts 
  mu_chain <- sigma2_chain <- matrix(0, nrow = nsamples, ncol = C+1)
  z_chain <- z0
  sigma2_mu <- sdpriormu^2

  # Pass these as parameters or guess based on something smart
  mu_current <- c(0, 1, 2)
  sigma2_current <- c(1, 1, 1)

  # Generate Markov Chain
  if(verbose) cat("Running Markov Chain iterations...", "\n")
  for(t in seq_len(nsamples)){
      # Update mus
      stts <- tapply(as.vector(y), INDEX = as.vector(z_chain), 
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
      z_chain <- inner_gibbs_conditional(z_chain,
          cond_weights = w, R = mrfi@Rmat, theta = theta, ncycles = 1)
      # TODO: add some counter for the latenf field in each pixel
      if(verbose) cat("\r", t)
  }
  if(verbose) cat("\n", "Done!", "\n")

  dfmus <- as.data.frame(mu_chain)
  colnames(dfmus) <- 0:C
  dfmus <- cbind(t = 1:nsamples, dfmus, par = "mu")
  
  dfsigma <- as.data.frame(sigma2_chain)
  colnames(dfsigma) <- 0:C
  dfsigma <- cbind(t = 1:nsamples, dfsigma, par = "sigma2")

  dfpars <- rbind(dfmus, dfsigma)
  
  out <- list(dfpars = dfpars, ll = llapprox)

  class(out) <- "hmrfbayes_out"
  return(out)
}
