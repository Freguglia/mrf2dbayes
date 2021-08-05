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
  cat(glue("mrfbayes output with {max(x$df$t)} observations from the posterior distribution."), "\n")
}

#' @rdname mrfbayes_out
#' @importFrom dplyr group_by summarize n
#' @importFrom stats quantile sd
#' @importFrom rlang .data
#' @export
summary.mrfbayes_out <- function(object, burnin = 0.25, ...){
  df <- object$df
  if(burnin < 1) burnin <- burnin*max(df$t)
  df <- df[df$t>burnin,]
  if(!object$rj) {
    stts <- summarize(group_by(df, .data$position, .data$interaction),
                      q025 = quantile(.data$value, probs = 0.025),
                      mean = mean(.data$value),
                      q975 = quantile(.data$value, probs = 0.975),
                      sd = sd(.data$value))
  } else {
    df <- df[df$value != 0, ]
    trange <- range(df$t)
    stts <- summarize(group_by(df, .data$position, .data$interaction),
                      prob = n()/diff(trange),
                      q025 = quantile(.data$value, probs = 0.025),
                      mean = mean(.data$value),
                      q975 = quantile(.data$value, probs = 0.975),
                      sd = sd(.data$value))
  }
  return(stts)
}

#' @rdname mrfbayes_out
#' @importFrom ggplot2 ggplot aes geom_line geom_rect theme_bw facet_wrap geom_hline
#' @importFrom dplyr %>% mutate group_by arrange ungroup
#' @export
plot.mrfbayes_out <- function(x, burnin = 0,
                              thin = 50, ...){
  if(!x$rj){
    stts <- summary(x, burnin = burnin)
    tmax <- max(x$df$t)
    p <- ggplot(x$df) +
      geom_line(aes(x = .data$t, y = .data$value, color = .data$position)) +
      geom_rect(data = stts,
                  aes(xmin = 0, xmax = tmax,
                      ymin = .data$q025, ymax = .data$q975,
                      fill = .data$position), alpha = 0.1) +
      theme_bw() +
      facet_wrap(~.data$interaction)

    return(p)
  } else {
    df <- x$df
    df <- df %>%
      group_by(.data$position, .data$interaction) %>%
      mutate(dif = c(0, diff(.data$t) - 1)) %>%
      mutate(grp = cumsum(.data$dif)) %>%
      ungroup() %>%
      filter(t %% thin == 0)
    p <- ggplot(df, aes(x = .data$t, y = .data$value, color = .data$position)) +
      geom_line(aes(group = interaction(.data$grp, .data$position))) +
      theme_bw() +
      geom_hline(yintercept = 0) +
      facet_wrap(~.data$interaction)

    return(p)
  }
}
