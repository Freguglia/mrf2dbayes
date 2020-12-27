#' @name mrfbayes_out
#' @title Output class of mrfbayes
#'
#' @param x a `mrfbayes_out` object.
#' @param object a `mrfbayes_out` object.
#' @param burnin number of (first) observations to exclude from calculations and plots.
#' @param ... other arguments (not used by this method).
#'
#' @importFrom glue glue

#' @rdname mrfbayes_out
#' @export
print.mrfbayes_out <- function(x, ...){
  cat(glue("mrfbayes output with {max(x$df$t)} observations from the posterior distribution.\n"))
}

#' @rdname mrfbayes_out
#' @importFrom dplyr group_by summarize
#' @importFrom stats quantile sd
#' @importFrom rlang .data
#' @export
summary.mrfbayes_out <- function(object, burnin = 0.25, ...){
  df <- object$df
  extra_args <- list(...)
  if(!is.null(extra_args$burnin)){
    burnin <- extra_args$burnin
    if(burnin < 1) burnin <- burnin*max(df$t)
    df <- df[df$t>burnin,]
  }
  stts <- summarize(group_by(df, .data$name),
                    q025 = quantile(.data$value, probs = 0.025),
                    mean = mean(.data$value),
                    q975 = quantile(.data$value, probs = 0.975),
                    sd = sd(.data$value))
  return(stts)
}

#' @rdname mrfbayes_out
#' @importFrom ggplot2 ggplot aes geom_line
#' @export
plot.mrfbayes_out <- function(x, ...){
  p <- ggplot(x$df, aes(x = .data$t, y = .data$value, color = .data$name)) +
    geom_line()
  return(p)
}
