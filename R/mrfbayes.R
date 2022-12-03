#' @name mrfbayes
#' @title Metropolis-Hasting algorithm for Markov Random Fields on lattices
#'
#' @inheritParams llapprox
#' @param z The observed random field.
#' @param llapprox A `llapprox` object based on the maximal interaction
#' structure.
#' @param nsamples Number of MCMC samples.
#' @param init_theta Initial values of the MCMC algorithm. Set to "zero" to
#'   automatically create a vector equal to zero with appropriate length or
#'   "pl" to use the Maximum Pseudolikelihood estimator of `z` as the initial
#'   value.
#' @param sdprior Sample Deviation of the Normal distributions used as prior.
#' @param sdkernel Sample Deviation of the Normal distributions of the
#'   transition kernel.
#'
#' @return An object of class `mrfbayes_out`
#'
#' @importFrom mrf2d expand_array smr_stat smr_array
#' @export
mrfbayes <- function(z, llapprox,
                     nsamples = 1000, init_theta = "zero",
                     sdprior = 10, sdkernel = 0.05,
                     verbose = interactive()){
  start_time <- Sys.time()
  # Initialize meta-parameters and do some checks
  family <- llapprox@family
  mrfi <- llapprox@mrfi
  T_zobs <- smr_stat(z, mrfi, family)
  if(llapprox@pass_entire){
    z_arg <- z
  } else {
    z_arg <- T_zobs
  }
  fdim <- length(T_zobs)
  C <- length(unique(as.vector(z))) - 1
  if(length(setdiff(0:C, unique(as.vector(z)))) > 1){
    stop("'z' should be a matrix with values in 0,...,C for some integer C.")
  }
  if(!class(llapprox) == "llapprox"){
    stop("'llapprox' must be an object created with the llapprox() function.")
  }

  # Initialize values
  if(is.character(init_theta)){
    if(init_theta == "zero"){
      current_theta <- T_zobs*0
    } else if(init_theta == "pl") {
      current_theta <- mrf2d::fit_pl(z, mrfi, family)$theta
      current_theta <- mrf2d::smr_array(current_theta, family)
    }
  } else if(is.array(init_theta)){

  } else {
    current_theta <- init_theta
  }

  resmat <- matrix(0, nrow = nsamples, ncol = fdim)

  # Run MCMC
  for(i in 1:nsamples){
    # Propose
    proposed_theta <- current_theta + rnorm(fdim, mean = 0, sd = sdkernel)

    # Compute acceptance probability
    logA <- llapprox@lafn(z_arg, proposed_theta) +
      sum(dnorm(proposed_theta, sd = sdprior, log = TRUE)) -
      llapprox@lafn(z_arg, current_theta) -
      sum(dnorm(current_theta, sd = sdprior, log = TRUE))
    if(log(runif(1)) < logA){
      current_theta <- proposed_theta
    }

    # Store results
    resmat[i,] <- current_theta

    if(verbose){
      cat("\r", i)
    }
  }
  if(verbose){
    cat("\n")
  }

  resdf <- as.data.frame(resmat)
  resdf$t <- 1:nsamples

  resdf <- tidyr::pivot_longer(resdf, cols = -"t")
  resdf <- cbind(resdf, mrf2d::vec_description(mrfi, family, C))
  resdf <- resdf[,c("t", "position", "interaction", "value")]

  end_time <- Sys.time()

  out <- list(df = tibble::as_tibble(resdf), ll = llapprox, rj = FALSE,
              ptime = as.numeric(difftime(end_time, start_time, units = "secs")))
  class(out) <- "mrfbayes_out"
  return(out)
}
