#' Bayesian calibration with a second predictor
#'
#' Fits a Bayesian regression of D47 on temperature and a second predictor (\code{Ion}) using Stan.
#'
#' @param calibrationData The calibration dataset. It needs the columns \code{Temperature} (10^6/T^2, with T in K), \code{TempError} (1 SD of \code{Temperature}), \code{D47}, \code{D47error} (1 SD of \code{D47}), and \code{Ion}.
#' @param IonError Ion measurement error (1 SD). Use a single value for all samples or one value per sample. Use \code{NULL} to treat \code{Ion} as known without error.
#' @param numSavedSteps Total number of posterior draws to keep across all chains.
#' @param chains Number of MCMC chains.
#' @param cores Number of cores used to run the chains in parallel.
#' @param MC Deprecated. Use \code{cores} instead.
#'
#' @details The model accounts for errors in \code{Temperature} and \code{D47}, and in \code{Ion} when \code{IonError} is given. Its \code{sigma} is the scatter left after removing the measurement error in \code{D47}. The true values of each predictor have a normal prior with a mean and SD estimated from the data. The intercept (alpha) and both slopes (beta and gamma) have flat priors.
#'
#' @return An object of class \code{bayclump_cal} with a single \code{stanfit} named \code{BLM1_fit}. Use \code{cal.draws()} to get the posterior draws.
#'
#' @import rstan
#'
#' @export

cal.ion.bayesian <- function(calibrationData,
                             IonError = NULL,
                             numSavedSteps = 3000,
                             chains = 4,
                             cores = getOption("mc.cores", 1L),
                             MC = NULL) {

  .check_cal_data(calibrationData, c("Temperature", "TempError", "D47", "D47error", "Ion"))

  # The MC argument is kept for older code
  if (!is.null(MC)) {
    warning("MC is deprecated. Use cores instead.")
    if (isTRUE(MC)) cores <- min(chains, parallel::detectCores(), na.rm = TRUE)
  }

  data <- list(
    N = nrow(calibrationData),
    x = calibrationData$Temperature,
    tau = calibrationData$TempError,
    ion = calibrationData$Ion,
    y = calibrationData$D47,
    y_err = calibrationData$D47error
  )

  # Start the true values at their measured values
  init <- function() list(x_true = data$x, mu_x = mean(data$x), sigma_x = sd(data$x))

  if (is.null(IonError)) {
    name <- "ion_noerrors"
  } else {
    # Use one error per sample, or the same error for all samples
    if (!length(IonError) %in% c(1, nrow(calibrationData)) || !all(IonError > 0)) {
      stop("IonError must be positive, with length 1 or one value per observation.")
    }
    data$ion_tau <- rep(IonError, length.out = nrow(calibrationData))
    name <- "ion_errors"
    init <- function() list(x_true = data$x, mu_x = mean(data$x), sigma_x = sd(data$x), ion_true = data$ion, mu_ion = mean(data$ion), sigma_ion = sd(data$ion))
  }

  fit <- .sample(name, data, c("alpha", "beta", "gamma", "sigma", "log_lik"), numSavedSteps, chains, cores, init = init)
  .finish_fit(list(BLM1_fit = fit), NULL, chains)
}
