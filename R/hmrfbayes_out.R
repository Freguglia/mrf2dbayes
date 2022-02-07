#' @rdname hmrfbayes
#' @param what Either `"pars"`, `"theta"` or `"zprobs"`.
#' @importFrom ggplot2 geom_tile scale_fill_gradient
#' @export
plot.hmrfbayes_out <- function(x, what = "pars", burnin = 0.25, ...){
  C <- x$ll@C
  if(what == "pars"){
    df <- tidyr::pivot_longer(x$dfpars, cols = as.character(0:C),
                              names_to = "k")
    p <- ggplot(df, aes(x = t, y = value, color = k)) +
      geom_line() +
      facet_wrap(~par, scales = "free") +
      theme_bw()
  } else if(what == "theta"){
		p <- ggplot(x$dftheta, aes(x = t, y = value, color = position)) +
      geom_line() +
      facet_wrap(~interaction) +
      theme_bw()
  } else if(what == "zprobs"){
    l <- lapply(seq(dim(x$zprobs)[3]), function(b) x$zprobs[ , , b])
    l <- lapply(l, function(slice){
                  data.frame(x = as.vector(row(slice)),
                             y = as.vector(col(slice)),
                             value = as.vector(slice))
          })
    df <- dplyr::bind_rows(l, .id = "k")
    df$k <- as.numeric(df$k) - 1
    p <- ggplot(df, aes(x = x, y = y, fill = value)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "black") +
      facet_wrap(~k) +
      theme_bw()
  } else {
    stop("'what' must be either 'pars', 'theta' or 'field'")
  }
	return(p)
}
