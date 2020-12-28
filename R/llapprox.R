#' @name llapprox-class
#' @title Log-Likelihood function approximations
#'
#' @description Objects for Log-Likelihood function approximations.
#'
#' @slot lafn The log-likelihood function approximation `lafn(z, theta)`.
#'   It is either a function of the complete lattice configuration or the
#'   vector of sufficient statistics, depending on the value of `pass_entire`.
#' @slot family a parameter restriction family from `mrf2d`.
#' @slot mrfi a `mrfi` object with the interaction structure.
#' @slot C Maximum integer value of the random fields considered
#'   (i.e., mrf support ranging from 0 to C)
#' @slot pass_entire `logical` value indicating whether the first argument of
#'   `lafn()` should be the entire lattice (`TRUE`) or the vector of sufficient
#'   statistics (`FALSE`).
#' @slot method The method used to create that function.
#' @slot internal_data Internal data to be used by the approximating function.
#' @exportClass llapprox
setClass("llapprox",
    representation(lafn = "function",
        family = "character",
        mrfi = "mrfi",
        C = "numeric",
        pass_entire = "logical",
        method = "character",
        internal_data = "list"))

setMethod("show", "llapprox",
    function(object){
        method_name <- switch(object@method,
            pseudo = "Pseudolikelihood",
            gt = "Geyer-Thompson method")
        cat(glue::glue("Log-Likelihood funcion approximation via {method_name}"), "\n")
    })

#' @rdname llapprox-class
#'
#'
#' @param refz A reference field to generate the approximation function. It is
#'   used to infer the lattice dimensions, support and any other relevant
#'   properties.
#' @param mrfi `mrfi` object.
#' @param family A parameter restriction family from `mrf2d`.
#' @param method Which type of approximation should be created. Options are:
#'   * `"pseudo"` for Pseudolikelihood approximation (which will result in
#'   Pseudoposteriors).
#'   * `"gt"` for Geyer-Thompson.
#'   * `"wanglandau"` for the Wang-Landau algorithm.
#' @param ... Extra options passed depending on the method used.
#'   * If the selected method is `"gt"`:
#'     * `reference`: A reference parameter value (theta), either as an array
#'     or a vector of appropriate length.
#'     * `nsamples`: The number of samples to be used in the approximation.
#'     * `ncycles`: The number of Gibbs Sampler cycles between each sample.
#'   * If the selected methid is `"wanglandau"`:
#'     * `reference`:
#' @param verbose `logical` value indicating wheter the algorithm progress
#'   should be printed.
#'
#'
#'
#' @return a `llapprox` object.
#'
#' @importFrom Brobdingnag as.brob
#'
#' @author Victor Freguglia
#' @export
llapprox <- function(refz, mrfi, family, method = "pseudo",
    verbose = interactive(), ...){
    la <- methods::new("llapprox")
    C <- length(unique(as.vector(refz))) - 1

    la@family <- family
    la@mrfi <- mrfi
    la@C <- C

    extra_args <- list(...)

    if(method == "pseudo"){
        la@method <- "pseudo"
        la@pass_entire <- TRUE
        la@lafn <- function(z, theta_vec){
            theta_arr <- mrf2d::expand_array(theta_vec, family, mrfi, C)
            mrf2d::pl_mrf2d(z, mrfi, theta_arr, log_scale = TRUE)
        }
    } else if(method == "gt"){
        la@method <- "gt"
        la@pass_entire <- FALSE

        if(is.null(extra_args$nsamples)){
            stop("The argument 'nsamples' must be specified for method 'gt'.")
        }
        if(is.null(extra_args$ncycles)){
            stop("The argument 'ncycles' must be specified for method 'gt'.")
        }

        # Get which theta to sample from
        if(is.null(extra_args$reference)){
            warning("'reference' is not specified in 'extra_args'. Using Maximum Pseudolikelihood estimator of 'refz' as reference value.")
            theta_ref_arr <- mrf2d::fit_pl(refz, mrfi, family)$theta
        } else {
            if(!is.array(extra_args$reference)){
                theta_ref_arr <- mrf2d::expand_array(extra_args$reference, family, mrfi, C)
            } else {
                theta_ref_arr <- extra_args$reference
            }
        }

        # Get sample size and number of cycles
        zmc <- mrf2d::rmrf2d(dim(refz), mrfi, theta_ref_arr, cycles = 60)
        Tzmc <- mrf2d::smr_stat(zmc, mrfi, family)
        if(verbose) cat("Sampling ergodic chain for Monte-Carlo approximations: \n")
        Tzmat <- mrf2d::rmrf2d_mc(zmc, mrfi, theta_ref_arr, family,
            nmc = extra_args$nsamples,
            burnin = 60, cycles = extra_args$ncycles)
        if(verbose) cat("\n")

        # Define log-likelihood function approximation
        theta_ref <- mrf2d::smr_array(theta_ref_arr, family)
        la@lafn <- function(z, theta_vec){
            Hs <- as.brob(as.vector(Tzmat %*% (theta_vec - theta_ref)))
            logzeta <- log(Brobdingnag::sum(exp(Hs))/length(Hs))
            return(as.vector(t(z) %*% theta_vec) - logzeta)
        }
        la@internal_data <- list(mcsamples = Tzmat)
    } else if(method == "wl") {
        # TODO: Implement Wang-Landau approximation
    } else {
        stop(glue::glue("'{method}' is not a valid method."))
    }

    return(la)
}
