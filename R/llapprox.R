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
#' @slot ptime Time (in seconds) used to compute the approximation.
#'
#' @importFrom glue glue
#' @exportClass llapprox
setClass("llapprox",
    representation(lafn = "function",
        family = "character",
        mrfi = "mrfi",
        C = "numeric",
        pass_entire = "logical",
        method = "character",
        internal_data = "list",
        ptime = "numeric"))

setMethod("show", "llapprox",
    function(object){
        method_name <- switch(object@method,
            pseudo = "Pseudolikelihood",
            gt = "Geyer-Thompson method")
        cat(glue("Log-Likelihood funcion approximation via {method_name}"), "\n")
        cat(glue("Interaction structure: {mrf2d::mrfi_to_string(object@mrfi)}"), "\n")
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
#'   * `"adj.pseudo"` for mode and curvature adjusted Pseudolikelihood
#'   approximation.
#'   * `"gt"` for Geyer-Thompson.
#'   * `"wanglandau"` for the Wang-Landau algorithm.
#' @param ... Extra options passed depending on the method used.
#'   * If the selected method is `"adj.pseudo"`:
#'     * `gamma_seq`: a sequence to be used by the Stochastic Approximation
#'     alforithm to find the Maximum Likelihood Estimator. See `?mrf2d::fit_sa`
#'     * `nsamples`: Number of samples drawn for the Monte-Carlo estimate of 
#'     the Likelihood function Hessian. Set to 1000 if unspecified.
#'   * If the selected method is `"gt"`:
#'     * `reference`: A reference parameter value (theta), either as an array
#'     or a vector of appropriate length.
#'     * `nsamples`: The number of samples to be used in the approximation.
#'     * `ncycles`: The number of Gibbs Sampler cycles between each sample.
#'   * If the selected method is `"wanglandau"`:
#'     * `theta_list`: A list of vectors with reference points for the 
#'     parameter vector.
#'     * `cycles`: Number of cycles in between samples. Defaults to 2.
#'     * `h`: Bandwidth to be used in the smoother function. Defaults to 10.
#' @param verbose `logical` value indicating wheter the algorithm progress
#'   should be printed.
#'
#'
#'
#' @return a `llapprox` object.
#'
#' @importFrom Brobdingnag as.brob cbrob
#' @importFrom mrf2d expand_array fit_pl pl_mrf2d smr_array rmrf2d
#' @importFrom purrr map2 map map_dbl map_chr
#'
#' @author Victor Freguglia
#' @export
llapprox <- function(refz, mrfi, family, method = "pseudo",
    verbose = interactive(), ...){
    start_time <- Sys.time()
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
            theta_arr <- expand_array(theta_vec, family, mrfi, C)
            pl_mrf2d(z, mrfi, theta_arr, log_scale = TRUE)
        }
    } else if(method == "adj.pseudo"){
        la@method <- "adj.pseudo"
        la@pass_entire <- TRUE

        if(is.null(extra_args$gamma_seq)){
            stop("The argument 'gamma_seq' must be specified for method 'adj.pseudo'.")
        }
        if(is.null(extra_args$nsamples)){
          nsamples <- 1000
        } else {
          nsamples <- extra_args$nsamples
        }
        
        # Find MLE and MPLE
        if(verbose) {cat("Obtaining Maximum Pseudolikelihood estimates...")}
        plfit <- fit_pl(refz, mrfi, family, return_optim = TRUE)
        mple <- smr_array(plfit$theta, family)
        if(verbose) {cat("Done!\n")}
        
        if(verbose) {cat("Obtaining Maximum Likelihood estimates via Stochastic Approximation...\n")}
        mlfit <- mrf2d::fit_sa(refz, mrfi, family, extra_args$gamma_seq,
                               verbose = verbose, init = mple)
        mle <- mrf2d::smr_array(mlfit$theta, family)
        if(verbose) {cat("Done!\n")}

        if(verbose) {cat("Generating Monte-Carlo samples to estimate Likelihood Hessian...\n")}
        samples <- as.matrix(mrf2d::rmrf2d_mc(refz, mrfi, mlfit$theta, family, nmc = nsamples,
                                    verbose = verbose))
        Tbar <- as.matrix(apply(samples, MARGIN = 2, mean))
        Etz <- apply(samples, MARGIN = 1, function(x) x%*%t(x)) %>% 
            as.matrix() %>% apply(MARGIN = 1, mean) %>% as.matrix()
        if(length(Tbar) == 1) Etz <- mean(Etz)
        Hsa <- -(Etz - (Tbar %*% t(Tbar)))
        Hpl <- plfit$opt.hessian
        N <- chol(as.matrix(-Hsa))
        M <- chol(as.matrix(-Hpl))
        W <- as.matrix(solve(M)%*%N)
        if(verbose) {cat("Done!\n")}

        la@lafn <- function(z, theta_vec){
            theta_arr <- mrf2d::expand_array(mple + as.vector(W%*%(theta_vec - mle)), family, mrfi, C)
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
            burnin = 60, cycles = extra_args$ncycles, verbose = verbose)
        # Define log-likelihood function approximation
        theta_ref <- mrf2d::smr_array(theta_ref_arr, family)
        la@lafn <- function(z, theta_vec){
            Hs <- as.brob(as.vector(Tzmat %*% (theta_vec - theta_ref)))
            logzeta <- log(Brobdingnag::sum(exp(Hs))/length(Hs))
            return(as.vector(t(z) %*% theta_vec) - logzeta)
        }
        la@internal_data <- list(mcsamples = Tzmat)

    } else if(method == "wl") {
        if(is.null(extra_args$theta_list)){
            stop("The argument 'theta_list' must be specified for method 'wl'.")
        }
        if(is.null(extra_args$gamma_seq)){
            stop("The argument 'gamma_seq' must be specified for method 'wl'.")
        }
        if(is.null(extra_args$cycles)){
          cycles <- 2
        } else {cycles <- extra_args$cycles}
        if(is.null(extra_args$h)){
          h <- 10
        } else {h <- extra_args$h}
        theta_list <- extra_args$theta_list
        gamma_seq <- extra_args$gamma_seq
        la@method <- "wl"
        la@pass_entire <- FALSE
        theta_is_array <- all(map_chr(theta_list, class) == "array")
        if(theta_is_array){
            theta_list_a <- theta_list
            theta_list_v <- map(theta_list, ~mrf2d::smr_array(.x, family))
        } else {
            theta_list_a <- map(theta_list, ~mrf2d::expand_array(.x, family = family, mrfi = mrfi, C = C))
            theta_list_v <- theta_list
        }
        kmix <- length(theta_list)
        N <- length(gamma_seq)
        theta_list_m <- do.call(rbind, theta_list_v)

        I <- integer(N)
        z <- rmrf2d(dim(refz), mrfi = mrfi, theta = theta_list_a[[sample(1:kmix, 1)]], cycles = 60)
        ss <- mrf2d::smr_stat(z, mrfi, family)
        ssmat <- matrix(nrow = N, ncol = length(ss))
        
        cvec <- numeric(kmix) + max(theta_list_m %*% ss)
        cvec_hist <- matrix(nrow = N, ncol = length(cvec)) 
        
        
        for(i in 1:N){
          I[i] <- sample(1:kmix, 1, 
                         prob = exp((theta_list_m %*% ss) - cvec))
          z <- rmrf2d(z, mrfi = mrfi, theta = theta_list_a[[I[i]]], cycles = cycles)
          ss <- smr_stat(z, mrfi, family)
          ssmat[i, ] <- ss
          cvec_hist[i, ] <- cvec - cvec[1]
          cvec <- cvec + gamma_seq[i]*((seq_len(kmix)==I[i]) - 1/kmix)
          cvec <- cvec + max((theta_list_m %*% ss) - cvec)
          if(verbose) cat("\r", i)
        }
        la@internal_data <- list(stat_hist = ssmat, cvec_hist = cvec_hist, I_hist = I)
        mixdf <- map(seq_along(theta_list), ~ssmat[I==.x,])
        
        smoother <- function(theta){
          dists <- map_dbl(theta_list_v, ~sum((.x - theta)^2))
          s <- exp(-0.5/h*dists)
          return(s/sum(s))
        }
        cvec2 <- colMeans(as.data.frame(cvec_hist[(N*0.75):N,]))
        zeta <- function(theta){
          Hs <- map2(mixdf, theta_list_v, ~as.matrix(.x) %*% (theta - .y))
          Hs <- map(Hs, ~exp(as.brob(as.vector(.x))))
          Hs <- map(Hs, ~Brobdingnag::sum(.x)/length(.x))
          Hs <- do.call(cbrob, Hs)
          Hs <- Hs * exp(as.brob(cvec2))
          weights <- smoother(theta)
          Hs <- Hs*weights
          Hs %>% Brobdingnag::sum() %>% log() %>% as.numeric()
        }

        la@lafn <- function(z, theta_vec){
            sum(z*theta_vec) - zeta(theta_vec) 
        }
        
    } else {
        stop(glue::glue("'{method}' is not a valid method."))
    }
    end_time <- Sys.time()
    la@ptime <- as.numeric(difftime(end_time, start_time, units = "secs"))
    return(la)
}
