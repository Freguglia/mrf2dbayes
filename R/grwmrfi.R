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
#' @export
grwmrfi <- function(z, llapprox, maximal_mrfi, init_mrfi = maximal_mrfi,
                    nsamples = 10000, nbatch = 100){

}
