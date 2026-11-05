#' Generate a dataset reflecting the priors used to run the analyses
#'
#' @param prior Either \code{Informative} or \code{Weak}, with the same values used in \code{cal.bayesian()}. The \code{Uninformative} option is a flat (improper) prior and can't be sampled.
#' @param n number of observations to simulate
#'
#' @importFrom stats rnorm
#'
#' @export


rec.prior <- function(prior, n = 1000) {

  if (!prior %in% c("Informative", "Weak")) {
    stop("Use `Informative` or `Weak`. `Uninformative` is a flat (improper) prior and can't be sampled.")
  }

  if (prior == "Informative") {
    params <- cbind.data.frame(
      parameter = c("alpha", "beta"),
      mean = c(0.231, 0.039),
      sd = c(0.065, 0.004)
    )
    params
  } else {
    params <- cbind.data.frame(
      parameter = c("alpha", "beta"),
      mean = c(0.231, 0.039),
      sd = c(1, 1)
    )
    params
  }
  x <- rnorm(n, 11, 0.0001)
  xC <- sqrt(10 ^ 6 / x) - 273.15

  data <- cbind.data.frame(alpha = rnorm(n, params[1, 2], params[1, 3]),
                           beta = rnorm(n, params[2, 2], params[2, 3]),
                           x,
                           xC)
  attr(data, "priors") <- prior
  attr(data, "params") <- params
  return(data)
}
