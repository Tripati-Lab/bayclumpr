#' Fit OLS regression models on a given calibration dataset
#'
#' Fits an ordinary least squares regression of D47 on temperature to bootstrap replicates of the calibration dataset. Measurement errors are not used.
#'
#' @param data The calibration dataset. It needs the columns \code{Temperature} (10^6/T^2, with T in K) and \code{D47}.
#' @param replicates Number of bootstrap replicates
#' @param samples Number of samples per bootstrap replicate. Defaults to the number of rows in \code{data}.
#'
#' @return A \code{data.frame} with one row per replicate and the columns \code{alpha} (intercept), \code{beta} (slope), and \code{sigma} (residual SD of D47 around the fitted line, including measurement error). Failed replicates are removed with a warning.
#'
#' @examples
#' ds <- cal.dataset(error = "S1", nobs = 50, seed = 1)
#' ols <- cal.ols(ds, replicates = 100)
#' colMeans(ols)
#'
#' @importFrom stats lm.fit
#'
#' @export

cal.ols <- function(data, replicates = 1000, samples = NULL) {
  .check_cal_data(data, c("Temperature", "D47"))
  .bootstrap(data, replicates, samples, function(d) {
    co <- unname(lm.fit(cbind(1, d$Temperature), d$D47)$coefficients)
    c(alpha = co[1], beta = co[2], sigma = .resid_sd(d, co[1], co[2]))
  })
}
