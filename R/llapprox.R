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
                                  pseudo = "Pseudolikelihood")
            cat(glue::glue("Log-Likelihood approximation via {method_name}\n"))
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
#' @param extra_args A list of extra options passed depending on the method
#'   used.
#'
#' @return a `llapprox` object.
#'
#' @author Victor Freguglia
#' @export
llapprox <- function(refz, mrfi, family, method = "pseudo", extra_args = NULL){
  la <- methods::new("llapprox")
  C <- length(unique(as.vector(refz))) - 1

  la@family <- family
  la@mrfi <- mrfi
  la@C <- C

  if(method == "pseudo"){
    la@method <- "pseudo"
    la@pass_entire <- TRUE
    la@lafn <- function(z, theta_vec){
      theta_arr <- mrf2d::expand_array(theta_vec, family, mrfi, C)
      mrf2d::pl_mrf2d(z, mrfi, theta_arr, log_scale = TRUE)
    }
  } else {
    stop(glue::glue("'{method}' is not a valid method."))
  }

  return(la)
}
