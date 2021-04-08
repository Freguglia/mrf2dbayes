#' @name grwmrfi
#' @title Graph Random Walk-based Metropolis-Hastings on Interaction Structures
#'
#' @description Metropolis-Hastings algorithm over the space of interaction
#' structures to be used as a Bayesian model selection tool.
#'
#' @param llapprox A `llapprox` object based on the maximal interaction
#' structure.
#' @param z The observed random field.
#' @param nsamples Number of MCMC samples.
#' @param nbatch Number of samples from the prior distribution in each
#' update of the acceptance ratios.
#' @param sdprior Sample Deviation of the Normal distributions used as prior.
#'
#' @return An object of class `grw_out`
#' @importFrom mrf2d mrfi_to_string
#' @export
grwmrfi <- function(z, llapprox, nsamples = 10000,
                    nbatch = 10, first_batch = 1000,
                    sdprior = 2, init_mrfi = mrfi(1),
                    verbose = interactive()){
  family <- llapprox@family
  Rmax <- llapprox@mrfi
  npos <- length(Rmax)
  T_zobs <- smr_stat(z, Rmax, family)
  if(llapprox@pass_entire){
    z_arg <- z
  } else {
    z_arg <- T_zobs
  }
  fdim <- length(T_zobs)
  dim_per_group <- fdim/npos

  current_mrfi <- init_mrfi
  mrfi_list <- mrfi_to_string(current_mrfi)
  included <- rep(TRUE, npos)
  current_mrfi_string <- mrfi_to_string(current_mrfi)
  psi <- ns <- mrfi_chain <- numeric(nsamples)
  psi <- as.brob(psi)

  for(t in seq_len(nsamples)){
    change <- sample(seq_len(npos),size = 1)
    change_pos <- Rmax[change]
    is_new <- FALSE
    if(length(current_mrfi + change_pos) > length(current_mrfi)){
         proposed_mrfi <- current_mrfi + change_pos
    } else { proposed_mrfi <- current_mrfi - change_pos }
    string_proposed <- mrfi_to_string(proposed_mrfi)
    if(!(string_proposed %in% mrfi_list)){
      mrfi_list <- c(mrfi_list, string_proposed)
      is_new <- TRUE
    }

    # Update psi and ns for current mrfi
    numtheta <- ifelse(is_new, first_batch, nbatch)
    which_current <- which(mrfi_list == current_mrfi_string)
    reduced_thetas <- lapply(seq_len(numtheta), function(x) rnorm(n = dim_per_group*npos, sd = sdprior)*rep(included, each = dim_per_group))
    current_liks <- unlist(lapply(reduced_thetas, function(x) llapprox@lafn(z_arg, x)))
    current_liks <- as.brob(current_liks)
    current_liks <- exp(current_liks)
    psi[which_current] <- (ns[which_current]*psi[which_current] + Brobdingnag::sum(current_liks))/ (ns[which_current] + nbatch)
    ns[which_current] <- ns[which_current] + numtheta

    # Update psi and ns for proposed mrfi
    which_proposed <- which(mrfi_list == string_proposed)
    included_proposed <- included
    included_proposed[change] <- !included_proposed[change]
    reduced_thetas <- lapply(seq_len(numtheta), function(x) rnorm(n = dim_per_group*npos, sd = sdprior)*rep(included_proposed, each = dim_per_group))
    current_liks <- unlist(lapply(reduced_thetas, function(x) llapprox@lafn(z_arg, x)))
    current_liks <- as.brob(current_liks)
    current_liks <- exp(current_liks)
    psi[which_proposed] <- (ns[which_proposed]*psi[which_proposed] + Brobdingnag::sum(current_liks))/ (ns[which_proposed] + nbatch)
    ns[which_proposed] <- ns[which_proposed] + numtheta

    U <- runif(1)
    if(U < psi[which_proposed]/(2^length(proposed_mrfi))/psi[which_current]*(2^length(current_mrfi))){
      mrfi_chain[t] <- which_proposed
      current_mrfi <- proposed_mrfi
      current_mrfi_string <- string_proposed
      included[change] <- !included[change]
    } else {
      mrfi_chain[t] <- which_current
    }

    if(verbose) cat("\r", t)
  }

  return(list(mrfi_list = mrfi_list, mrfi_chain = mrfi_chain, ns = ns, psi = psi))
}
