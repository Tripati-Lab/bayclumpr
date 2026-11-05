#' Fit York regression models on a given calibration dataset
#'
#' Fits a York regression of D47 on temperature to bootstrap replicates of the calibration dataset, using \code{IsoplotR::york()}. The errors in both variables are used for each observation.
#'
#' @param data The calibration dataset. It needs the columns \code{Temperature} (10^6/T^2, with T in K), \code{TempError} (1 SD of \code{Temperature}), \code{D47}, and \code{D47error} (1 SD of \code{D47}).
#' @param replicates Number of bootstrap replicates
#' @param samples Number of samples per bootstrap replicate. Defaults to the number of rows in \code{data}.
#'
#' @return A \code{data.frame} with one row per replicate and the columns \code{alpha} (intercept), \code{beta} (slope), and \code{sigma} (residual SD of D47 around the fitted line, including measurement error). Failed replicates are removed with a warning.
#'
#' @examples
#' ds <- cal.dataset(error = "S1", nobs = 50, seed = 1)
#' yrk <- cal.york(ds, replicates = 100)
#' colMeans(yrk)
#'
#' @importFrom IsoplotR york
#'
#' @export

cal.york <- function(data, replicates = 1000, samples = NULL) {
  .check_cal_data(data, c("Temperature", "TempError", "D47", "D47error"))
  .bootstrap(data, replicates, samples, function(d) {
    fit <- york(cbind(d$Temperature, d$TempError, d$D47, d$D47error))
    a <- unname(fit$a[1])
    b <- unname(fit$b[1])
    c(alpha = a, beta = b, sigma = .resid_sd(d, a, b))
  })
}
