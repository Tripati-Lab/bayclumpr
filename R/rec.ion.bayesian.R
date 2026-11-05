#' Bayesian temperature reconstructions with a second predictor
#'
#' Reconstructs temperatures (in degrees C) from the posterior draws of \code{cal.ion.bayesian()}, the measured D47, and the \code{Ion} value of each sample.
#'
#' @param calModel Posterior draws from \code{cal.draws()}, a result from \code{cal.ion.bayesian()}, or the \code{stanfit} in it.
#' @param recData The reconstruction dataset. It needs the columns \code{Sample}, \code{D47}, \code{D47error} (1 SD of the D47 value), and \code{Ion}. Add \code{N} (number of replicates) or \code{D47sd} (SD of the replicates) to also get \code{SD}. The output of \code{avg.replicates()} has both. It also needs \code{IonError} (1 SD of \code{Ion}) when \code{useIonError = TRUE}.
#' @param ndraws Number of draws used for each sample.
#' @param useIonError Whether to account for the error in \code{Ion}.
#' @param level Probability covered by the interval.
#' @param ... Arguments from earlier versions (\code{iter}, \code{postcalsamples}, and \code{MC}). They are ignored.
#'
#' @details With flat priors on 10^6/T^2 and \code{Ion}, the posterior of \code{Ion} follows its measurement, and the posterior of 10^6/T^2 for each draw is normal. Both are sampled directly, so no MCMC run is needed. Draws with 10^6/T^2 at or below zero are not physically possible. They are removed with a warning.
#'
#' The uncertainty includes the calibration parameters, the residual SD of the calibration (\code{sigma}), \code{D47error}, and \code{IonError} when used.
#'
#' @return A \code{data.frame} with one row per sample. \code{meanTemp} is the posterior mean temperature (in degrees C). \code{SE} is the standard error of the temperature, propagated from \code{D47error}. \code{SD} is the SD of the temperature of a single replicate, propagated from the SD of the replicates (\code{D47error} x sqrt(\code{N}), or \code{D47sd} when \code{N} is missing). It is \code{NA} without replicate information. Both include the uncertainty in the calibration. \code{median}, \code{lower}, and \code{upper} are the posterior median and the limits of the interval, based on \code{SE}.
#'
#' @examples
#' \donttest{
#' ds <- cal.dataset(error = "S1", nobs = 50, seed = 1)
#' ds$Ion <- rnorm(nrow(ds))
#' fit <- cal.ion.bayesian(ds, IonError = 0.1, numSavedSteps = 1000)
#' recData <- data.frame(Sample = "A", D47 = 0.6, D47error = 0.01, Ion = 0, IonError = 0.1)
#' rec.ion.bayesian(fit, recData, useIonError = TRUE)
#' }
#'
#' @importFrom stats rnorm sd
#'
#' @export

rec.ion.bayesian <- function(calModel, recData, ndraws = 10000, useIonError = FALSE, level = 0.95, ...) {
  old <- names(list(...))
  if (length(old) > 0) {
    warning("These arguments are no longer used and were ignored: ", paste(old, collapse = ", "), ".")
  }

  draws <- if (is.data.frame(calModel)) calModel else cal.draws(calModel)
  if (!"gamma" %in% names(draws)) {
    stop("The calibration needs a slope for Ion (gamma). Use a model from cal.ion.bayesian().")
  }

  cols <- c("Sample", "D47", "D47error", "Ion")
  if (useIonError) cols <- c(cols, "IonError")
  .check_rec_data(recData, cols)
  res <- .rec_exact(draws, recData, ndraws, level = level, useIon = TRUE, useIonError = useIonError)
  attr(res, "plug") <- NULL
  res
}
