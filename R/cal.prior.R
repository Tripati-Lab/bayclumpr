#' Generate a dataset reflecting the priors used to run the analyses
#'
#' @param prior The priors on alpha and beta. Use \code{"Informative"} or \code{"Weak"}, a list such as \code{list(alpha = c(0.231, 0.065), beta = c(0.039, 0.004))}, or a result from \code{cal.bayesian()} to use its priors. Flat priors are improper and can't be sampled.
#' @param n number of observations to simulate
#'
#' @return A \code{data.frame} with \code{n} draws of \code{alpha} and \code{beta}.
#'
#' @importFrom stats rnorm
#'
#' @export

cal.prior <- function(prior, n = 1000) {
  params <- if (inherits(prior, "bayclump_cal")) attr(prior, "priors") else .get_priors(prior)
  if (!is.list(params)) {
    stop("Flat priors are improper and can't be sampled.")
  }

  data <- cbind.data.frame(
    alpha = rnorm(n, params$alpha[1], params$alpha[2]),
    beta = rnorm(n, params$beta[1], params$beta[2])
  )
  attr(data, "priors") <- if (is.character(prior)) prior else "Custom"
  attr(data, "params") <- cbind.data.frame(
    parameter = c("alpha", "beta"),
    mean = c(params$alpha[1], params$beta[1]),
    sd = c(params$alpha[2], params$beta[2])
  )
  data
}
