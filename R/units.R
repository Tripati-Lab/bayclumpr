#' Convert temperatures to and from 10^6/T^2
#'
#' \code{temp2x()} converts temperatures in degrees C or kelvin to 10^6/T^2 (T in K), the scale used by all calibrations in \code{bayclumpr}. \code{x2temp()} converts them back.
#'
#' @param temp Temperatures in degrees C or kelvin.
#' @param tempError Optional 1 SD of \code{temp}, in the same units.
#' @param x Temperatures as 10^6/T^2 (T in K).
#' @param xError Optional 1 SD of \code{x}.
#' @param unit Unit of \code{temp}, or the unit returned by \code{x2temp()}: \code{"C"} or \code{"K"}.
#'
#' @details Errors are propagated to first order. The SD of 10^6/T^2 is 2 x 10^6 x SD(T) / T^3, with T in kelvin. This approximation works well when the error is small relative to T, which is the case for calibration data.
#'
#' @return \code{temp2x()} returns a \code{data.frame} with \code{Temperature} and, when \code{tempError} is given, \code{TempError}. These names match the columns needed by the calibration functions. \code{x2temp()} returns a \code{data.frame} with \code{Temp} and, when \code{xError} is given, \code{TempSD}.
#'
#' @examples
#' temp2x(c(5, 25, 50), tempError = c(0.5, 0.5, 1))
#' x2temp(c(12.9, 11.3, 9.6), xError = 0.02)
#'
#' @export

temp2x <- function(temp, tempError = NULL, unit = c("C", "K")) {
  unit <- match.arg(unit)
  tempK <- if (unit == "C") temp + 273.15 else temp
  if (any(tempK <= 0)) {
    stop("Temperatures must be above absolute zero.")
  }

  out <- data.frame(Temperature = 10^6 / tempK^2)
  if (!is.null(tempError)) out$TempError <- 2 * 10^6 * tempError / tempK^3
  out
}

#' @rdname temp2x
#' @export

x2temp <- function(x, xError = NULL, unit = c("C", "K")) {
  unit <- match.arg(unit)
  if (any(x <= 0)) {
    stop("x must be positive.")
  }

  tempK <- sqrt(10^6 / x)
  out <- data.frame(Temp = if (unit == "C") tempK - 273.15 else tempK)
  if (!is.null(xError)) out$TempSD <- xError * tempK^3 / (2 * 10^6)
  out
}
