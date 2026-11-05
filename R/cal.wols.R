#' Fit weighted OLS regression models on a given calibration dataset
#'
#' Fits a weighted least squares regression of D47 on temperature to bootstrap replicates of the calibration dataset. The weights don't come from the reported errors. They come from a regression of the absolute OLS residuals on the fitted values, so observations in regions with more scatter get less weight.
#'
#' @param data The calibration dataset. It needs the columns \code{Temperature} (10^6/T^2, with T in K) and \code{D47}.
#' @param replicates Number of bootstrap replicates
#' @param samples Number of samples per bootstrap replicate. Defaults to the number of rows in \code{data}.
#'
#' @details A replicate fails when the regression of the residuals predicts a non-positive value, since the weights can't be computed.
#'
#' @return A \code{data.frame} with one row per replicate and the columns \code{alpha} (intercept), \code{beta} (slope), and \code{sigma} (residual SD of D47 around the fitted line, including measurement error). Failed replicates are removed with a warning.
#'
#' @examples
#' ds <- cal.dataset(error = "S1", nobs = 50, seed = 1)
#' wols <- cal.wols(ds, replicates = 100)
#' colMeans(wols)
#'
#' @importFrom stats lm.fit lm.wfit
#'
#' @export

cal.wols <- function(data, replicates = 1000, samples = NULL) {
  .check_cal_data(data, c("Temperature", "D47"))
  .bootstrap(data, replicates, samples, function(d) {
    X <- cbind(1, d$Temperature)
    ols <- lm.fit(X, d$D47)

    # Model the size of the residuals to get the weights
    s <- lm.fit(cbind(1, d$D47 - ols$residuals), abs(ols$residuals))$fitted.values
    if (any(s <= 0)) stop("Non-positive weights.")

    co <- unname(lm.wfit(X, d$D47, w = 1 / s^2)$coefficients)
    c(alpha = co[1], beta = co[2], sigma = .resid_sd(d, co[1], co[2]))
  })
}
