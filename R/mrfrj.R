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
#' @importFrom dplyr filter
#' @export
mrfrj <- function(z, llapprox,
                  nsamples = 1000, init_theta = "zero",
                  sdprior = 1, sdkernel = 0.005,
                  sdjump = 0.05, sdcenter = 0.3,
                  logpenalty = log(prod(dim(z))),
                  verbose = interactive()){
  start_time <- Sys.time()

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
      current_theta <- rnorm(length(T_zobs), mean = 0, sd = 0.01)*0
    } else if(init_theta == "pl") {
      if(!is.null(llapprox@internal_data$mple)){
        current_theta <- llapprox@internal_data$mple
      } else {
        current_theta <- mrf2d::fit_pl(z, maximal_mrfi, family)$theta
        current_theta <- mrf2d::smr_array(current_theta, family)
      }
    }
  } else if(is.array(init_theta)){
    current_theta <- init_theta
    current_theta <- mrf2d::smr_array(current_theta, family)
  } else {
    current_theta <- init_theta
  }
  current_lafn <- llapprox@lafn(z_arg, current_theta)

  resmat <- matrix(0, nrow = nsamples, ncol = fdim)
  included <- rep(FALSE, length(maximal_mrfi))
  proposals <- data.frame(t = 1:nsamples,
                          move = factor(character(nsamples),
                                        levels = c("swap", "within", "jump", "share", "center")),
                          logA = numeric(nsamples),
                          extraInfo = character(nsamples))

  # Run MCMC
  for(i in 1:nsamples){
    # Propose move
    move <- sample(c("jump", "within", "share", "swap", "center"), size = 1)
    proposals$move[i] <- move

    # Walk within the current model
    if(move == "within"){
      proposed_theta <- current_theta + rnorm(fdim, mean = 0, sd = sdkernel)*(current_theta!=0)
  proposed_lafn <- llapprox@lafn(z_arg, proposed_theta)

      # Compute acceptance probability
      logA <- proposed_lafn +
        sum(dnorm(proposed_theta, sd = sdprior, log = TRUE)) -
        current_lafn -
        sum(dnorm(current_theta, sd = sdprior, log = TRUE))
      proposals$logA[i] <- logA

      if(log(runif(1)) < logA){
        current_theta <- proposed_theta
        current_lafn <- proposed_lafn
      }

    # Add or delete a position
    } else if(move == "jump"){
      jump_pos <- sample(seq_len(npos), size = 1)
      vec_jump <- seq_len(npos) == jump_pos

      if(included[jump_pos]){ # Delete the selected position
        proposed_theta <- current_theta * rep(included*(!vec_jump), each = dim_per_group)
        proposed_lafn <- llapprox@lafn(z_arg, proposed_theta)

        logA <- proposed_lafn +
          sum(dnorm(proposed_theta, sd = sdprior, log = TRUE)) -
          current_lafn -
          sum(dnorm(current_theta, sd = sdprior, log = TRUE)) +
          sum(dnorm(current_theta[rep(vec_jump, each = dim_per_group)], sd = sdjump, log = TRUE)) +
          logpenalty
        proposals$logA[i] <- logA
        proposals$extraInfo[i] <- "delete"
        if(log(runif(1)) < logA){
          current_theta <- proposed_theta
          current_lafn <- proposed_lafn
          included <- included - vec_jump
        }

      } else { # Add the selected position
        proposed_theta <- current_theta + rnorm(fdim, mean = 0, sd = sdjump)*rep(vec_jump, each = dim_per_group)
        proposed_lafn <- llapprox@lafn(z_arg, proposed_theta)

        logA <- proposed_lafn +
          sum(dnorm(proposed_theta, sd = sdprior, log = TRUE)) -
          current_lafn -
          sum(dnorm(current_theta, sd = sdprior, log = TRUE)) -
          sum(dnorm(proposed_theta[rep(vec_jump, each = dim_per_group)], sd = sdjump, log = TRUE)) -
          logpenalty
        proposals$logA[i] <- logA
        proposals$extraInfo[i] <- "add"

        if(log(runif(1)) < logA){
          current_theta <- proposed_theta
          current_lafn <- proposed_lafn
          included <- included + vec_jump
        }
      }

    # Share weights of a particular interaction between two positions
    } else if(move == "share"){
      if(sum(included) > 1){
        # Propose
        to_share <- sample(which(as.logical(included)), 2, replace = FALSE)
        w <- runif(1, 0, 1)
        norm1 <- sqrt(sum(current_theta[((to_share[1] - 1)*dim_per_group + 1):((to_share[1])*dim_per_group)]^2))
        norm2 <- sqrt(sum(current_theta[((to_share[2] - 1)*dim_per_group + 1):((to_share[2])*dim_per_group)]^2))
        sum_norms <- norm1 + norm2
        proposed_theta <- current_theta
        proposed_theta[((to_share[1] - 1)*dim_per_group + 1):((to_share[1])*dim_per_group)] <-
          proposed_theta[((to_share[1] - 1)*dim_per_group + 1):((to_share[1])*dim_per_group)]/norm1*w*sum_norms
        proposed_theta[((to_share[2] - 1)*dim_per_group + 1):((to_share[2])*dim_per_group)] <-
        proposed_theta[((to_share[2] - 1)*dim_per_group + 1):((to_share[2])*dim_per_group)]/norm2*(1-2)*sum_norms

        proposed_lafn <- llapprox@lafn(z_arg, proposed_theta)

        logA <- proposed_lafn +
          sum(dnorm(proposed_theta, sd = sdprior, log = TRUE)) -
          current_lafn -
          sum(dnorm(current_theta, sd = sdprior, log = TRUE))
        proposals$logA[i] <- logA
        if(is.nan(logA)) logA <- -Inf
        if(log(runif(1)) < logA){
          current_theta <- proposed_theta
          current_lafn <- proposed_lafn
        }
      }
    } else if(move == "center"){
      if(sum(included) > 1){
        to_center <- sample(which(as.logical(included)), 2, replace = FALSE)
        proposed_theta <- current_theta
        idx1 <- ((to_center[1] - 1)*dim_per_group + 1):((to_center[1])*dim_per_group)
        idx2 <- ((to_center[2] - 1)*dim_per_group + 1):((to_center[2])*dim_per_group)
        idx_int <- sample(1:dim_per_group)
        idx1 <- idx1[idx_int]
        idx2 <- idx2[idx_int]

        u <- rnorm(1, mean = 0, sd = sdcenter)
        proposed_theta[idx1] <- proposed_theta[idx1] + u
        proposed_theta[idx2] <- proposed_theta[idx2] - u
        proposed_lafn <- llapprox@lafn(z_arg, proposed_theta)

        logA <- proposed_lafn +
          sum(dnorm(proposed_theta, sd = sdprior, log = TRUE)) -
          current_lafn -
          sum(dnorm(current_theta, sd = sdprior, log = TRUE))
        proposals$logA[i] <- logA
        if(is.nan(logA)) logA <- -Inf
        if(log(runif(1)) < logA){
          current_theta <- proposed_theta
          current_lafn <- proposed_lafn
        }
      }
    } else if(move == "swap"){
      if(sum(included) > 1 && sum(included) < length(included)){
        to_remove <- sample(which(as.logical(included)), 1)
        to_include <- sample(which(!as.logical(included)), 1)
        proposed_theta <- current_theta

        proposed_theta[((to_include - 1)*dim_per_group + 1):(to_include*dim_per_group)] <-
          proposed_theta[((to_remove - 1)*dim_per_group + 1):(to_remove*dim_per_group)]
        proposed_theta[((to_remove - 1)*dim_per_group + 1):(to_remove*dim_per_group)] <- 0

        proposed_lafn <- llapprox@lafn(z_arg, proposed_theta)

        logA <- proposed_lafn -
          current_lafn
        proposals$logA[i] <- logA
        if(is.nan(logA)) logA <- -Inf
        if(log(runif(1)) < logA){
          current_theta <- proposed_theta
          current_lafn <- proposed_lafn
          included[to_remove] <- FALSE
          included[to_include] <- TRUE
        }
      }
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

  resdf <- tidyr::gather(resdf, "key", "value", -t) %>%
    arrange(t)
  resdf <- cbind(resdf, mrf2d::vec_description(maximal_mrfi, family, C))
  resdf <- resdf[,c("t", "position", "interaction", "value")]
  resdf <- resdf[resdf$value != 0.0,]
  resdf <- tibble::as_tibble(resdf)

  end_time <- Sys.time()
  ptime <- as.numeric(difftime(end_time, start_time, units = "secs"))

  out <- list(df = resdf, ll = llapprox, rj = TRUE,
              ptime = ptime, proposals = proposals)
  class(out) <- "mrfbayes_out"
  return(out)
}
