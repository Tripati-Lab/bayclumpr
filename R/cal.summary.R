#' Summarize a calibration
#'
#' Summarizes the intercept, slope, and residual SD of a calibration. It works the same way for bootstrap calibrations (\code{cal.ols()} and others) and Bayesian calibrations (\code{cal.bayesian()} and \code{cal.ion.bayesian()}).
#'
#' @param fit A \code{data.frame} of bootstrap replicates or posterior draws, a result from \code{cal.bayesian()} or \code{cal.ion.bayesian()}, or a single \code{stanfit} from it.
#' @param model The model to use when \code{fit} holds more than one, such as \code{"errors"}.
#' @param level Probability covered by the interval.
#'
#' @details \code{SE} is the standard error of each estimate. It is the SD of the bootstrap replicates or of the posterior draws. \code{sigma} is itself an SD: the scatter of individual samples around the calibration line. Its \code{SE} describes how well that scatter is known.
#'
#' For Bayesian fits, \code{ESS} is the effective sample size (\code{n_eff}) and \code{Rhat} the potential scale reduction factor of each parameter, as reported by \code{rstan}. Look for \code{Rhat} at or below 1.01 and at least 100 effective draws per chain. They are \code{NA} for bootstrap replicates, which are independent draws, and for a \code{data.frame} of posterior draws, which no longer records the chains.
#'
#' @return A \code{data.frame} with one row per parameter and the columns \code{parameter}, \code{estimate} (mean), \code{SE}, \code{median}, \code{lower}, \code{upper}, \code{ESS}, and \code{Rhat}. Mixed models also have a \code{Material} column.
#'
#' @examples
#' ds <- cal.dataset(error = "S1", nobs = 50, seed = 1)
#' ols <- cal.ols(ds, replicates = 200)
#' cal.summary(ols)
#'
#' @importFrom stats quantile sd
#'
#' @export

cal.summary <- function(fit, model = NULL, level = 0.95) {
  draws <- if (is.data.frame(fit)) fit else cal.draws(fit, model)
  pars <- intersect(c("alpha", "beta", "gamma", "sigma"), names(draws))
  probs <- c(0.5, (1 - level) / 2, 1 - (1 - level) / 2)

  # ESS and Rhat need the chains, so they come from the stanfit
  diag <- NULL
  if (!is.data.frame(fit)) {
    sf <- .get_stanfit(fit, model)
    diag <- summary(sf, pars = intersect(pars, sf@sim$pars_oi))$summary[, c("n_eff", "Rhat"), drop = FALSE]
  }

  # One row per parameter. In mixed models, alpha and beta are indexed by material and sigma is shared.
  summ <- function(d, j = NULL) {
    do.call(rbind, lapply(pars, function(p) {
      q <- quantile(d[[p]], probs, names = FALSE)
      key <- if (is.null(j) || p == "sigma") p else paste0(p, "[", j, "]")
      ess <- if (is.null(diag)) NA else diag[key, "n_eff"]
      rhat <- if (is.null(diag)) NA else diag[key, "Rhat"]
      data.frame(parameter = p, estimate = mean(d[[p]]), SE = sd(d[[p]]), median = q[1], lower = q[2], upper = q[3], ESS = ess, Rhat = rhat)
    }))
  }

  # Mixed models are summarized by material, in the order used by the fit
  if ("Material" %in% names(draws)) {
    mats <- unique(draws$Material)
    out <- do.call(rbind, lapply(seq_along(mats), function(j) cbind(Material = mats[j], summ(draws[draws$Material == mats[j], ], j))))
  } else {
    out <- summ(draws)
  }
  row.names(out) <- NULL
  out
}
