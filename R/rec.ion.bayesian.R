#' Bayesian temperature reconstructions with a second predictor
#'
#' Reconstructs temperatures (in degrees C) from the posterior draws of \code{cal.ion.bayesian()}, the measured D47, and the \code{Ion} value of each sample.
#'
#' @param calModel Posterior draws from \code{cal.draws()}, a result from \code{cal.ion.bayesian()}, or the \code{stanfit} in it.
#' @param recData The reconstruction dataset. It needs the columns \code{Sample}, \code{D47}, \code{D47error} (1 SD of the D47 value), and \code{Ion}. It also needs \code{IonError} (1 SD of \code{Ion}) when \code{useIonError = TRUE}.
#' @param ndraws Number of draws used for each sample.
#' @param useIonError Whether to account for the error in \code{Ion}.
#' @param ... Arguments from earlier versions (\code{iter}, \code{postcalsamples}, and \code{MC}). They are ignored.
#'
#' @details With flat priors on 10^6/T^2 and \code{Ion}, the posterior of \code{Ion} follows its measurement, and the posterior of 10^6/T^2 for each draw is normal. Both are sampled directly, so no MCMC run is needed. Draws with 10^6/T^2 at or below zero are not physically possible. They are removed with a warning.
#'
#' @return A \code{data.frame} with \code{meanTemp} (posterior mean, in degrees C) and \code{error} (1 SD) for each sample. The error includes uncertainty in the calibration parameters, the residual SD of the calibration (\code{sigma}), \code{D47error}, and \code{IonError} when used.
#'
#' @importFrom stats rnorm sd
#'
#' @export

rec.ion.bayesian <- function(calModel, recData, ndraws = 10000, useIonError = FALSE, ...) {
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
  .rec_exact(draws, recData, ndraws, useIon = TRUE, useIonError = useIonError)
}
