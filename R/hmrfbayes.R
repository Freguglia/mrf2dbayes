#' @export
hmrfbayes <- function(y, llapprox, nsamples = 1000, z0 = NULL, sdpriormu = 10, alpha1 = 0.001, alpha2 = 0.001){
  C <- llapprox@C
  mrfi <- llapprox@mrfi

  # Check and set initial field config
  if(is.null(z0)){
    z0 <- matrix(sample(0:2, size = prod(dim(y)), replace = TRUE),
                 nrow = nrow(y), ncol = ncol(y))
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
          dnorm(y, mean = mu_current[x], sd = sqrt(sigma2_current[x]), log = TRUE)
      })
      w <- array(do.call(c, w2), dim = c(dim(y), length(mu_current)))
      z_chain <- inner_gibbs_conditional(z_chain,
          cond_weights = w, R = mrfi@Rmat, theta = theta, ncycles = 1)
      # TODO: add some counter for the latenf field in each pixel
      cat("\r", t)
  }

  return(list(mu_chain, sigma2_chain))

}
