#' Bayesian temperature reconstructions
#'
#' Reconstructs temperatures (in degrees C) from the posterior draws of a Bayesian calibration and the measured D47 of each sample. This approach is congruent with the one used in McClelland et al. (2022).
#'
#' @param calModel Posterior draws from \code{cal.draws()}, a result from \code{cal.bayesian()}, or a single \code{stanfit} from it.
#' @param recData The reconstruction dataset. It needs the columns \code{Sample}, \code{D47}, and \code{D47error} (1 SD of the D47 value). If D47 is the mean of several replicates, \code{D47error} should be its standard error. Mixed models also need \code{Material}, with the same labels used in the calibration.
#' @param ndraws Number of draws used for each sample.
#' @param model The model to use when \code{calModel} holds more than one, such as \code{"errors"}.
#' @param ... Arguments from earlier versions (\code{iter}, \code{mixed}, \code{postcalsamples}, and \code{MC}). They are ignored.
#'
#' @details With a flat prior on 10^6/T^2, the posterior of 10^6/T^2 for each calibration draw is normal. Its mean is (D47 - alpha) / beta and its SD is sqrt(sigma^2 + D47error^2) / beta. This function samples it directly, so no MCMC run is needed. Draws with 10^6/T^2 at or below zero are not physically possible. They are removed with a warning.
#'
#' @return A \code{data.frame} with \code{meanTemp} (posterior mean, in degrees C) and \code{error} (1 SD) for each sample. The error includes uncertainty in the calibration parameters, the residual SD of the calibration (\code{sigma}), and \code{D47error}. For the errors model, \code{sigma} excludes the measurement error of the calibration samples. For the other two models, it includes it.
#'
#' @importFrom stats rnorm sd
#'
#' @export

rec.bayesian <- function(calModel, recData, ndraws = 10000, model = NULL, ...) {
  old <- names(list(...))
  if (length(old) > 0) {
    warning("These arguments are no longer used and were ignored: ", paste(old, collapse = ", "), ".")
  }

  draws <- if (is.data.frame(calModel)) calModel else cal.draws(calModel, model)
  .check_rec_data(recData, c("Sample", "D47", "D47error"))
  .rec_exact(draws, recData, ndraws)
}
