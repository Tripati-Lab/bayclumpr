#' Calibration line with its interval
#'
#' Computes the median calibration line and its interval over a range of temperatures. This is mostly used for plotting.
#'
#' @param data A \code{data.frame} with the columns \code{alpha} and \code{beta}, such as a result from \code{cal.ols()} or \code{cal.draws()}. For mixed models, use the rows of a single material.
#' @param from The lowest temperature (10^6/T^2).
#' @param to The highest temperature (10^6/T^2).
#' @param length.out Number of temperatures between \code{from} and \code{to}.
#' @param level Probability covered by the interval.
#'
#' @return A \code{data.frame} with the columns \code{x} (temperature), \code{median_est} (median D47), and \code{ci_lower_est} and \code{ci_upper_est} (limits of the interval).
#'
#' @examples
#' ds <- cal.dataset(error = "S1", nobs = 50, seed = 1)
#' ols <- cal.ols(ds, replicates = 100)
#' ci <- cal.ci(ols, from = min(ds$Temperature), to = max(ds$Temperature))
#' head(ci)
#'
#' @importFrom stats quantile
#'
#' @export

cal.ci <- function(data, from, to, length.out = 100, level = 0.95) {
  if ("Material" %in% names(data) && length(unique(data$Material)) > 1) {
    stop("data has several materials. Use the rows of a single material.")
  }

  # One row per temperature and one column per draw
  x <- seq(from, to, length.out = length.out)
  pred <- outer(x, data$beta) + matrix(data$alpha, nrow = length(x), ncol = nrow(data), byrow = TRUE)

  probs <- c(0.5, (1 - level) / 2, 1 - (1 - level) / 2)
  q <- t(apply(pred, 1, quantile, probs = probs, names = FALSE))

  data.frame(x = x, median_est = q[, 1], ci_lower_est = q[, 2], ci_upper_est = q[, 3])
}
