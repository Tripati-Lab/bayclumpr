#' Summarize a calibration
#'
#' Summarizes the intercept, slope, and residual SD of a calibration. It works the same way for bootstrap calibrations (\code{cal.ols()} and others) and Bayesian calibrations (\code{cal.bayesian()} and \code{cal.ion.bayesian()}).
#'
#' @param fit A \code{data.frame} of bootstrap replicates or posterior draws, a result from \code{cal.bayesian()} or \code{cal.ion.bayesian()}, or a single \code{stanfit} from it.
#' @param model The model to use when \code{fit} holds more than one, such as \code{"errors"}.
#' @param level Probability covered by the interval.
#'
#' @details \code{SE} is the standard error of each estimate. It is the SD of the bootstrap replicates or of the posterior draws. \code{sigma} is itself an SD: the scatter of individual samples around the calibration line. Its \code{SE} describes how well that scatter is known.
#'
#' @return A \code{data.frame} with one row per parameter and the columns \code{parameter}, \code{estimate} (mean), \code{SE}, \code{median}, \code{lower}, and \code{upper}. Mixed models also have a \code{Material} column.
#'
#' @examples
#' ds <- cal.dataset(error = "S1", nobs = 50, seed = 1)
#' ols <- cal.ols(ds, replicates = 200)
#' cal.summary(ols)
#'
#' @importFrom stats quantile sd
#'
#' @export

cal.summary <- function(fit, model = NULL, level = 0.95) {
  draws <- if (is.data.frame(fit)) fit else cal.draws(fit, model)
  pars <- intersect(c("alpha", "beta", "gamma", "sigma"), names(draws))
  probs <- c(0.5, (1 - level) / 2, 1 - (1 - level) / 2)

  # One row per parameter
  summ <- function(d) {
    do.call(rbind, lapply(pars, function(p) {
      q <- quantile(d[[p]], probs, names = FALSE)
      data.frame(parameter = p, estimate = mean(d[[p]]), SE = sd(d[[p]]), median = q[1], lower = q[2], upper = q[3])
    }))
  }

  # Mixed models are summarized by material
  if ("Material" %in% names(draws)) {
    mats <- unique(draws$Material)
    out <- do.call(rbind, lapply(mats, function(m) cbind(Material = m, summ(draws[draws$Material == m, ]))))
  } else {
    out <- summ(draws)
  }
  row.names(out) <- NULL
  out
}
