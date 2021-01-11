#' @rdname hmrfbayes_out
#' @param what Either `"pars"`, `"theta"` or `"field"`.
#' @importFrom ggplot2 ggplot aes geom_line geom_rect theme_bw facet_wrap
#' @export
plot.hmrfbayes_out <- function(x, burnin = 0.25, what = "pars", ...){
  C <- x$ll@C
  if(what == "pars"){
    df <- tidyr::pivot_longer(x$dfpars, cols = as.character(0:C),
                              names_to = "k")
    p <- ggplot(df, aes(x = t, y = value, color = k)) +
      geom_line() + 
      facet_wrap(~par, scales = "free")
    return(p)
  } else if(what == "theta"){

  } else if(what == "field"){

  } else {
    stop("'what' must be either 'pars', 'theta' or 'field'")
  }
}
