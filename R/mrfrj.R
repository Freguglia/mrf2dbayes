#' @name mrfrj
#' @title Reversible Jump for MRFs over different interaction structures
#'
#' @description A Metropolis-Hastings algorithm that allows jumps over spaces
#' of different dimensions. This algorithm allows spaces with restrictions
#' (zero equality) to be visitted, providing a model selection framework along
#' with the Bayesian solution for inference within models.
#'
#' @param llapprox A `llapprox` object based on the maximal interaction
#' structure.
#' @param z The observed random field.
#' @param llapprox The likelihood approximation to be used.
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
mrfrj <- function(z, llapprox,
                  nsamples = 1000, init_theta = "pl",
                  sdprior = 10, sdkernel = 0.01,
                  sdjump = 0.1,
                  logpenalty = log(prod(dim(z))),
                  verbose = interactive()){

  # Initialize meta-parameters and do some checks
  family <- llapprox@family
  maximal_mrfi <- llapprox@mrfi
  T_zobs <- smr_stat(z, maximal_mrfi, family)
  if(llapprox@pass_entire){
    z_arg <- z
  } else {
    z_arg <- T_zobs
  }
  fdim <- length(T_zobs)
  npos <- length(maximal_mrfi)
  dim_per_group <- fdim/npos
  C <- length(unique(as.vector(z))) - 1
  if(length(setdiff(0:C, unique(as.vector(z)))) > 1){
    stop("'z' should be a matrix with values in 0,...,C for some integer C.")
  }
  if(!class(llapprox) == "llapprox"){
    stop("'llapprox' must be an object created with the llapprox() function.")
  }

  # Initialize values
  ## Start in the maximal model
  if(is.character(init_theta)){
    if(init_theta == "zero"){
      current_theta <- T_zobs*0
    } else if(init_theta == "pl") {
      current_theta <- mrf2d::fit_pl(z, maximal_mrfi, family)$theta
      current_theta <- mrf2d::smr_array(current_theta, family)
    }
  } else if(is.array(init_theta)){

  } else {
    current_theta <- init_theta
  }

  resmat <- matrix(0, nrow = nsamples, ncol = fdim)
  included <- rep(TRUE, length(maximal_mrfi))

  # Run MCMC
  for(i in 1:nsamples){
    # Propose move
    move <- sample(c("jump", "within"), size = 1)

    # Walk within the current model
    if(move == "within"){
      proposed_theta <- current_theta + rnorm(fdim, mean = 0, sd = sdkernel)*(current_theta!=0)

      # Compute acceptance probability
      logA <- llapprox@lafn(z_arg, proposed_theta) +
        sum(dnorm(proposed_theta, sd = sdprior, log = TRUE)) -
        llapprox@lafn(z_arg, current_theta) -
        sum(dnorm(current_theta, sd = sdprior, log = TRUE))
      if(log(runif(1)) < logA){
        current_theta <- proposed_theta
      }

    # Add or delete a position
    } else if(move == "jump"){
      jump_pos <- sample(seq_len(npos), size = 1)
      vec_jump <- seq_len(npos) == jump_pos

      if(included[jump_pos]){ # Delete the selected position
        proposed_theta <- current_theta * rep(included*(!vec_jump), each = dim_per_group)

        logA <- llapprox@lafn(z_arg, proposed_theta) +
          sum(dnorm(proposed_theta, sd = sdprior, log = TRUE)) -
          llapprox@lafn(z_arg, current_theta) -
          sum(dnorm(current_theta, sd = sdprior, log = TRUE)) +
          sum(dnorm(current_theta[rep(vec_jump, each = dim_per_group)], sd = sdjump, log = TRUE)) +
          logpenalty
        if(log(runif(1)) < logA){
          current_theta <- proposed_theta
          included <- included - vec_jump
        }

      } else { # Add the selected position
        proposed_theta <- current_theta + rnorm(fdim, mean = 0, sd = sdjump)*rep(vec_jump, each = dim_per_group)

        logA <- llapprox@lafn(z_arg, proposed_theta) +
          sum(dnorm(proposed_theta, sd = sdprior, log = TRUE)) -
          llapprox@lafn(z_arg, current_theta) -
          sum(dnorm(current_theta, sd = sdprior, log = TRUE)) -
          sum(dnorm(proposed_theta[rep(vec_jump, each = dim_per_group)], sd = sdjump, log = TRUE)) -
          logpenalty

        if(log(runif(1)) < logA){
          current_theta <- proposed_theta
          included <- included + vec_jump
        }
      }

    # Shift a position
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
  resdf <- cbind(resdf, mrf2d::vec_description(maximal_mrfi, family, C))
  resdf <- resdf[,c("t", "position", "interaction", "value")]

  out <- list(df = tibble::as_tibble(resdf), ll = llapprox, rj = TRUE)
  class(out) <- "mrfbayes_out"
  return(out)
}
