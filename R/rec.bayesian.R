#' Bayesian temperature reconstructions
#'
#' Reconstructs temperatures (in degrees C) from the posterior draws of a Bayesian calibration and the measured D47 of each sample. This approach is congruent with the one used in McClelland et al. (2022).
#'
#' @param calModel Posterior draws from \code{cal.draws()}, a result from \code{cal.bayesian()}, or a single \code{stanfit} from it.
#' @param recData The reconstruction dataset. It needs the columns \code{Sample}, \code{D47}, and \code{D47error} (1 SD of the D47 value). If D47 is the mean of several replicates, \code{D47error} should be its standard error. Add \code{N} (number of replicates) or \code{D47sd} (SD of the replicates) to also get \code{SD}. The output of \code{avg.replicates()} has both. Mixed models also need \code{Material}, with the same labels used in the calibration.
#' @param ndraws Number of draws used for each sample.
#' @param model The model to use when \code{calModel} holds more than one, such as \code{"errors"}.
#' @param newMaterial For mixed models, whether the samples belong to a material that is not in the calibration. When \code{TRUE}, \code{Material} in \code{recData} is ignored. When \code{calModel} is a table of draws, get them with \code{cal.draws(fit, newMaterial = TRUE)} instead.
#' @param level Probability covered by the interval.
#' @param ... Arguments from earlier versions (\code{iter}, \code{mixed}, \code{postcalsamples}, and \code{MC}). They are ignored.
#'
#' @details With a flat prior on 10^6/T^2, the posterior of 10^6/T^2 for each calibration draw is normal. Its mean is (D47 - alpha) / beta and its SD is sqrt(sigma^2 + D47error^2) / beta. This function samples it directly, so no MCMC run is needed. Draws with 10^6/T^2 at or below zero are not physically possible. They are removed with a warning.
#'
#' The uncertainty includes the calibration parameters, the residual SD of the calibration (\code{sigma}), and \code{D47error}. For models with measurement errors, \code{sigma} excludes the measurement error of the calibration samples. For models without them, it includes it.
#'
#' @return A \code{data.frame} with one row per sample. \code{meanTemp} is the posterior mean temperature (in degrees C). \code{SE} is the standard error of the temperature, propagated from \code{D47error}. \code{SD} is the SD of the temperature of a single replicate, propagated from the SD of the replicates (\code{D47error} x sqrt(\code{N}), or \code{D47sd} when \code{N} is missing). It is \code{NA} without replicate information. Both include the uncertainty in the calibration. \code{median}, \code{lower}, and \code{upper} are the posterior median and the limits of the interval, based on \code{SE}.
#'
#' @examples
#' \donttest{
#' ds <- cal.dataset(error = "S1", nobs = 50, seed = 1)
#' fit <- cal.bayesian(ds, models = "errors", numSavedSteps = 1000)
#' recData <- data.frame(Sample = c("A", "B"), D47 = c(0.6, 0.7), D47error = 0.01, N = 3)
#' rec.bayesian(fit, recData)
#' }
#'
#' @importFrom stats rnorm sd
#'
#' @export

rec.bayesian <- function(calModel, recData, ndraws = 10000, model = NULL, newMaterial = FALSE, level = 0.95, ...) {
  old <- names(list(...))
  if (length(old) > 0) {
    warning("These arguments are no longer used and were ignored: ", paste(old, collapse = ", "), ".")
  }

  draws <- if (is.data.frame(calModel)) calModel else cal.draws(calModel, model, newMaterial = newMaterial)
  .check_rec_data(recData, c("Sample", "D47", "D47error"))
  res <- .rec_exact(draws, recData, ndraws, level = level)
  attr(res, "plug") <- NULL
  res
}
