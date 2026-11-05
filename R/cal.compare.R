#' Compare Bayesian calibration models
#'
#' Compares the models in a \code{cal.bayesian()} result using approximate leave-one-out cross-validation (PSIS-LOO).
#'
#' @param fit A result from \code{cal.bayesian()} with at least two models.
#'
#' @details In models with measurement errors, the pointwise log-likelihood is conditional on the latent true temperatures. A warning is given when any Pareto k value is above 0.7, since the comparison may then be unreliable.
#'
#' @return The result of \code{loo::loo_compare()}. The best model is in the first row.
#'
#' @examples
#' \donttest{
#' ds <- cal.dataset(error = "S1", nobs = 50, seed = 1)
#' fit <- cal.bayesian(ds, models = c("errors", "noerrors"), numSavedSteps = 1000)
#' cal.compare(fit)
#' }
#'
#' @importFrom loo extract_log_lik relative_eff loo_compare
#'
#' @export

cal.compare <- function(fit) {
  if (!inherits(fit, "bayclump_cal") || length(fit) < 2) {
    stop("cal.compare() needs a result from cal.bayesian() with at least two models.")
  }

  loos <- lapply(fit, function(f) {
    log_lik <- extract_log_lik(f, merge_chains = FALSE)
    loo::loo(log_lik, r_eff = relative_eff(log_lik))
  })

  # Flag unreliable PSIS-LOO estimates
  k <- unlist(lapply(loos, loo::pareto_k_values))
  if (any(k > 0.7)) {
    warning("Some Pareto k values are above 0.7, so the comparison may be unreliable.")
  }

  loo_compare(loos)
}
