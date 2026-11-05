#' Fit Deming regression models on a given calibration dataset
#'
#' Fits a Deming regression of D47 on temperature to bootstrap replicates of the calibration dataset, using \code{deming::deming()}. The errors in both variables are used for each observation.
#'
#' @param data The calibration dataset. It needs the columns \code{Temperature} (10^6/T^2, with T in K), \code{TempError} (1 SD of \code{Temperature}), \code{D47}, and \code{D47error} (1 SD of \code{D47}).
#' @param replicates Number of bootstrap replicates
#' @param samples Number of samples per bootstrap replicate. Defaults to the number of rows in \code{data}.
#'
#' @return A \code{data.frame} with one row per replicate and the columns \code{alpha} (intercept), \code{beta} (slope), and \code{sigma} (residual SD of D47 around the fitted line, including measurement error). Failed replicates are removed with a warning.
#'
#' @examples
#' ds <- cal.dataset(error = "S1", nobs = 50, seed = 1)
#' dem <- cal.deming(ds, replicates = 100)
#' colMeans(dem)
#'
#' @importFrom deming deming
#'
#' @export

cal.deming <- function(data, replicates = 1000, samples = NULL) {
  .check_cal_data(data, c("Temperature", "TempError", "D47", "D47error"))
  x_SE <- y_SE <- NULL
  .bootstrap(data, replicates, samples, function(d) {
    d$x_SE <- d$TempError
    d$y_SE <- d$D47error

    # The jackknife SEs are not needed, so they are skipped
    co <- unname(deming(D47 ~ Temperature, d, xstd = x_SE, ystd = y_SE, jackknife = FALSE)$coefficients)
    c(alpha = co[1], beta = co[2], sigma = .resid_sd(d, co[1], co[2]))
  })
}
