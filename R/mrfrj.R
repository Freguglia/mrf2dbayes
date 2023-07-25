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
#' @importFrom dplyr filter pull
#' @importFrom tidyr unite
#' @import data.table
#' @export
mrfrj <- function(z, llapprox,
                  nsamples = 1000, init_theta = "zero",
                  init_included = "zero",
                  sdprior = 1, sdkernel = 0.005,
                  sdbirth = 0.05,
                  kernel_probs = c(4,1,1,1,1),
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
  if(is.character(init_included)){
    if(init_included == "zero"){
      included <- rep(FALSE, length(maximal_mrfi))
    } else if (init_included == "full"){
      included <- rep(TRUE, length(maximal_mrfi))
    } else {
      stop("Invalid 'init_included'.")
    }
  } else {
    if(length(init_included) == length(maximal_mrfi)){
      included <- init_included
    } else {
      stop("'init_included' must have the same length as the maximal mrfi.")
    }
  }

  current_theta <- mrfrj_init_theta(init_theta, T_zobs, sdkernel,
                                    included, dim_per_group, llapprox,
                                    maximal_mrfi, family)
  current_lafn <- llapprox@lafn(z_arg, current_theta)

  resmat <- matrix(0, nrow = nsamples, ncol = fdim)

  move_list <- c("within",
                 "swap",
                 "death",
                 "birth",
                 "jump"
                 )

  proposals <- data.frame(t = 1:nsamples,
                          move = factor(character(nsamples),
                                        levels = move_list),
                          logA = numeric(nsamples),
                          logPP = numeric(nsamples),
                          extraInfo = character(nsamples))

  # Run MCMC
  for(i in 1:nsamples){
    # Propose move
    move <- sample(move_list, size = 1, prob = kernel_probs)
    proposals$move[i] <- move

    # Random walk proposal within the current model
    if(move == "within"){
      proposal <- mrfrj_step_within(current_theta, current_lafn, z_arg,
                                    fdim, sdkernel, sdprior, llapprox)
      proposals$logA[i] <- proposal$logA
      current_theta <- proposal$theta
      current_lafn <- proposal$lafn

    # Add or delete a position
    } else if(move == "birth"){
      if(sum(included) < npos){
        pos_to_add <- sample(which(!included), size = 1)
        theta_to_add <- rnorm(dim_per_group, mean = 0, sd = sdbirth)

        proposed_theta <- current_theta
        proposed_theta[((pos_to_add-1)*dim_per_group + 1):(pos_to_add*dim_per_group)] <- theta_to_add

        if(sum(included) > 0){
          wts <- rgamma(sum(included), shape = 1/10)
          wts <- wts/sum(wts)
          mult <- rep(included, each = dim_per_group)
          mult[mult] <- mult[mult]*rep(wts, each = dim_per_group)
          proposed_theta <- proposed_theta - theta_to_add*mult
        }

        proposed_lafn <- llapprox@lafn(z_arg, proposed_theta)

        logA <- proposed_lafn +
          sum(dnorm(proposed_theta[proposed_theta!=0], sd = sdprior, log = TRUE)) -
          current_lafn -
          sum(dnorm(current_theta[current_theta!=0], sd = sdprior, log = TRUE)) -
          logpenalty -
          sum(dnorm(theta_to_add, sd = sdbirth, log = TRUE)) +
          log(kernel_probs[3] * 1/(sum(included) + 1)) -
          log(kernel_probs[4] * 1/(npos-sum(included)))
        proposals$logA[i] <- logA

        if(log(runif(1)) < logA) {
          current_theta <- proposed_theta
          included[pos_to_add] <- TRUE
          current_lafn <- proposed_lafn
        }
      } else { # Could not add more positions, full model
        logA <- -Inf
        proposals$logA[i] <- logA
      }
    } else if(move == "death"){
      if(sum(included) > 0){
        pos_to_remove <- sample(which(included), size = 1)
        theta_to_remove <- current_theta[((pos_to_remove-1)*dim_per_group + 1):(pos_to_remove*dim_per_group)]
        remaining <- included - (seq_len(npos) == pos_to_remove)

        proposed_theta <- current_theta
        proposed_theta[((pos_to_remove-1)*dim_per_group + 1):(pos_to_remove*dim_per_group)] <- 0

        if(sum(remaining) > 0){
          wts <- rgamma(sum(remaining), shape = 1/10)
          wts <- wts/sum(wts)
          mult <- rep(remaining, each = dim_per_group)
          mult[mult] <- mult[mult]*rep(wts, each = dim_per_group)
          proposed_theta <- proposed_theta + theta_to_remove*mult
        }

        proposed_lafn <- llapprox@lafn(z_arg, proposed_theta)

        logA <- proposed_lafn +
          sum(dnorm(proposed_theta[proposed_theta!=0], sd = sdprior, log = TRUE)) -
          current_lafn -
          sum(dnorm(current_theta[current_theta!=0], sd = sdprior, log = TRUE)) +
          logpenalty +
          sum(dnorm(theta_to_remove, sd = sdbirth, log = TRUE)) +
          log(kernel_probs[4] * 1/(npos-sum(included)+1)) -
          log(kernel_probs[3] * 1/sum(included))
        proposals$logA[i] <- logA

        if(log(runif(1)) < logA) {
          current_theta <- proposed_theta
          included[pos_to_remove] <- FALSE
          current_lafn <- proposed_lafn
        }
      } else {
        logA <- -Inf
        proposals$logA[i] <- logA
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
      } else {
        logA <- -Inf
        proposals$logA[i] <- logA
      }
      # Add or delete a position
    } else if(move == "jump"){
      jump_pos <- sample(seq_len(npos), size = 1)
      vec_jump <- seq_len(npos) == jump_pos

      if(included[jump_pos]){ # Delete the selected position
        proposed_theta <- current_theta * rep(included*(!vec_jump), each = dim_per_group)
        proposed_lafn <- llapprox@lafn(z_arg, proposed_theta)

        logA <- proposed_lafn +
          sum(dnorm(proposed_theta[proposed_theta!=0], sd = sdprior, log = TRUE)) -
          current_lafn -
          sum(dnorm(current_theta[current_theta != 0], sd = sdprior, log = TRUE)) +
          sum(dnorm(current_theta[rep(vec_jump, each = dim_per_group)], sd = sdbirth, log = TRUE)) +
          logpenalty
        proposals$logA[i] <- logA
        proposals$extraInfo[i] <- "delete"
        if(log(runif(1)) < logA){
          current_theta <- proposed_theta
          current_lafn <- proposed_lafn
          included <- as.logical(included - vec_jump)
        }

      } else { # Add the selected position
        proposed_theta <- current_theta + rnorm(fdim, mean = 0, sd = sdbirth)*rep(vec_jump, each = dim_per_group)
        proposed_lafn <- llapprox@lafn(z_arg, proposed_theta)

        logA <- proposed_lafn +
          sum(dnorm(proposed_theta[proposed_theta!=0], sd = sdprior, log = TRUE)) -
          current_lafn -
          sum(dnorm(current_theta[proposed_theta!=0], sd = sdprior, log = TRUE)) -
          sum(dnorm(proposed_theta[rep(vec_jump, each = dim_per_group)], sd = sdbirth, log = TRUE)) -
          logpenalty
        proposals$logA[i] <- logA
        proposals$extraInfo[i] <- "add"

        if(log(runif(1)) < logA){
          current_theta <- proposed_theta
          current_lafn <- proposed_lafn
          included <- as.logical(included + vec_jump)
        }
      }
    }

    # Store results
    resmat[i,] <- current_theta
    proposals$logPP[i] <- current_lafn +
      sum(dnorm(current_theta[current_theta!=0], sd = sdprior, log = TRUE))

    if(verbose){
      cat("\r", i)
    }
  }
  if(verbose){
    cat("\n")
  }
  gc()

  colnames(resmat) <- mrf2d::vec_description(maximal_mrfi, family, C) %>%
    unite(name, position, interaction) %>%
    pull(name)

  resdf <- as.data.table(resmat)
  rm(resmat)
  resdf[,t:=1:nsamples]
  gc()
  resdf <- data.table::melt(as.data.table(resdf), id.vars = "t")[
    order(t)][
      value != 0][
        , c("position", "interaction") := tstrsplit(variable, "_", fixed=TRUE)][,variable:=NULL]
  gc()
  resdf <- tibble::as_tibble(resdf)

  end_time <- Sys.time()
  ptime <- as.numeric(difftime(end_time, start_time, units = "secs"))

  out <- list(df = resdf, ll = llapprox, rj = TRUE,
              ptime = ptime, proposals = proposals,
              last_included = included, last_theta = current_theta)
  class(out) <- "mrfbayes_out"
  return(out)
}

mrfrj_init_theta <- function(init_theta,
                             T_zobs,
                             sdkernel,
                             included,
                             dim_per_group,
                             llapprox,
                             maximal_mrfi,
                             family){
  if(is.character(init_theta)){
    if(init_theta == "zero"){
      current_theta <- rnorm(length(T_zobs), mean = 0, sd = sdkernel)*rep(included, each = dim_per_group)
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
}


mrfrj_step_within <- function(current_theta, current_lafn, z_arg, fdim, sdkernel, sdprior, llapprox){
  proposed_theta <- current_theta + rnorm(fdim, mean = 0, sd = sdkernel)*(current_theta!=0)
  proposed_lafn <- llapprox@lafn(z_arg, proposed_theta)

  logA <- proposed_lafn +
    sum(dnorm(proposed_theta[proposed_theta!=0], sd = sdprior, log = TRUE)) -
    current_lafn -
    sum(dnorm(current_theta[current_theta!=0], sd = sdprior, log = TRUE))

  if(log(runif(1)) < logA){
    next_theta <- proposed_theta
    next_lafn <- proposed_lafn
  } else {
    next_theta <- current_theta
    next_lafn <- current_lafn
  }
  return(list(theta = next_theta, lafn = next_lafn, logA = logA))
}
