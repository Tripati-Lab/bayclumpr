#' Temperature reconstructions from bootstrap calibrations
#'
#' Reconstructs temperatures (in degrees C) from the bootstrap replicates of a calibration on 10^6/T^2 (T in K) and the measured D47 of each sample.
#'
#' @param recData Reconstruction dataset. It needs the columns \code{Sample}, \code{D47}, and \code{D47error} (1 SD of the D47 value). If D47 is the mean of several replicates, \code{D47error} should be its standard error. Add \code{N} (number of replicates) or \code{D47sd} (SD of the replicates) to also get \code{SD}. The output of \code{avg.replicates()} has both.
#' @param obCal A \code{data.frame} with the slopes and intercepts of the calibration, such as a result from \code{cal.ols()}. It also needs \code{sigma} when \code{prediction = TRUE}.
#' @param ndraws Number of Monte Carlo draws used for each sample.
#' @param prediction Whether to include the residual scatter of the calibration (\code{sigma}). Use \code{TRUE} to get the same kind of uncertainty as \code{rec.bayesian()}.
#' @param level Probability covered by the interval.
#'
#' @details The uncertainty always includes \code{D47error} and the joint spread of slopes and intercepts in \code{obCal}. With \code{prediction = FALSE}, it describes the temperature predicted by the calibration line. With \code{prediction = TRUE}, it also includes the scatter of individual samples around that line. Draws with 10^6/T^2 at or below zero are not physically possible. They are removed with a warning.
#'
#' @return A \code{data.frame} with one row per sample. \code{meanTemp} is the temperature from the mean slope and intercept (in degrees C). \code{SE} is the standard error of the temperature, propagated from \code{D47error}. \code{SD} is the SD of the temperature of a single replicate, propagated from the SD of the replicates (\code{D47error} x sqrt(\code{N}), or \code{D47sd} when \code{N} is missing). It is \code{NA} without replicate information. Both include the uncertainty in the calibration. \code{median}, \code{lower}, and \code{upper} are the median and the limits of the interval, based on \code{SE}.
#'
#' @examples
#' ds <- cal.dataset(error = "S1", nobs = 50, seed = 1)
#' ols <- cal.ols(ds, replicates = 100)
#' recData <- data.frame(Sample = c("A", "B"), D47 = c(0.6, 0.7), D47error = 0.01, N = 3)
#' rec.clumped(recData, ols)
#' rec.clumped(recData, ols, prediction = TRUE)
#'
#' @export

rec.clumped <- function(recData,
                        obCal,
                        ndraws = 10000,
                        prediction = FALSE,
                        level = 0.95) {

  .check_rec_data(recData, c("Sample", "D47", "D47error"))

  if (prediction && is.null(obCal$sigma)) {
    stop("obCal needs a sigma column when prediction = TRUE. Refit the calibration with this version of bayclumpr.")
  }
  if (!prediction) obCal$sigma <- 0

  res <- .rec_exact(obCal, recData, ndraws, level = level)

  # Use the temperature from the mean slope and intercept as the point estimate
  res$meanTemp <- attr(res, "plug")
  attr(res, "plug") <- NULL
  res
}
