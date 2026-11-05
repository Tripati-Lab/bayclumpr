#' Generate a dataset reflecting the priors used to run the analyses
#'
#' Same as \code{cal.prior()}, with two extra columns: \code{x}, values of 10^6/T^2 close to 11, and \code{xC}, the same values in degrees C.
#'
#' @param prior The priors on alpha and beta. See \code{cal.prior()}.
#' @param n number of observations to simulate
#'
#' @return A \code{data.frame} with \code{n} rows and the columns \code{alpha}, \code{beta}, \code{x}, and \code{xC}.
#'
#' @examples
#' head(rec.prior("Weak", n = 100))
#'
#' @importFrom stats rnorm
#'
#' @export

rec.prior <- function(prior, n = 1000) {
  data <- cal.prior(prior, n)
  data$x <- rnorm(n, 11, 0.0001)
  data$xC <- sqrt(10 ^ 6 / data$x) - 273.15
  data
}
