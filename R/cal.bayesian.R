#' Bayesian calibration of the clumped isotope paleothermometer
#'
#' Fits up to four Bayesian linear regressions of D47 on temperature using Stan. The models cover every combination of two choices: whether measurement errors are used, and whether intercepts and slopes vary by material.
#'
#' @param calibrationData The calibration dataset. It needs the columns \code{Temperature} (10^6/T^2, with T in K), \code{TempError} (1 SD of \code{Temperature}), \code{D47}, and \code{D47error} (1 SD of \code{D47}). Mixed models also need \code{Material}.
#' @param models The models to fit. Use any of \code{"errors"}, \code{"noerrors"}, \code{"mixed_errors"}, and \code{"mixed_noerrors"}. \code{"mixed"} is kept as another name for \code{"mixed_noerrors"}.
#' @param priors Priors on the intercept (alpha) and slope (beta). Use \code{"Informative"}, \code{"Weak"}, or \code{"Flat"}, or a list such as \code{list(alpha = c(0.231, 0.065), beta = c(0.039, 0.004))} with the mean and SD of each normal prior. \code{"Uninformative"} is kept as another name for \code{"Flat"}.
#' @param numSavedSteps Total number of posterior draws to keep across all chains.
#' @param chains Number of MCMC chains.
#' @param cores Number of cores used to run the chains in parallel.
#' @param MC Deprecated. Use \code{cores} instead.
#' @param x A \code{bayclump_cal} object.
#' @param digits Number of significant digits to print.
#' @param ... Not used.
#'
#' @details The four models are:
#' \tabular{llll}{
#'   \code{models} \tab Name in the result \tab Measurement errors \tab Varies by \code{Material} \cr
#'   \code{"errors"} \tab \code{BLM1_fit} \tab Yes \tab No \cr
#'   \code{"noerrors"} \tab \code{BLM1_fit_NoErrors} \tab No \tab No \cr
#'   \code{"mixed_errors"} \tab \code{BLM3_fit_Errors} \tab Yes \tab Yes \cr
#'   \code{"mixed_noerrors"} \tab \code{BLM3_fit} \tab No \tab Yes
#' }
#'
#' Models with measurement errors use \code{TempError} and \code{D47error} for each observation. Their \code{sigma} is the scatter left after removing the measurement error in \code{D47}. The true temperatures have a normal prior with a mean and SD estimated from the data. Models without measurement errors ignore both errors, so their \code{sigma} includes the measurement error.
#'
#' Mixed models let intercepts and slopes vary by \code{Material} around shared means (partial pooling). The priors set by \code{priors} apply to the shared means. The SDs among materials have half-normal priors with scales of 0.1 for the intercept at the mean calibration temperature and 0.01 for the slope. With fewer than five materials, these SDs are mostly set by their priors. Mixed models are skipped when there are fewer than two materials.
#'
#' Temperature is centered inside every model to improve sampling. Results are reported for the usual intercept at 10^6/T^2 = 0, where the priors are also set.
#'
#' The \code{"Informative"} values are the ones used in earlier versions of this package. Check that they suit the reference frame of your data, or set your own.
#'
#' A warning is given when a model shows signs of poor convergence: Rhat above 1.01, fewer than 100 effective draws per chain, or divergent transitions.
#'
#' @return An object of class \code{bayclump_cal}. It is a list of \code{stanfit} objects named as in the table above, with the priors and a convergence summary stored as attributes. Printing it shows \code{cal.summary()} for each model and the convergence summary. Use \code{cal.summary()} to get these tables, \code{cal.draws()} to get the posterior draws, and \code{cal.compare()} to compare models.
#'
#' @examples
#' \donttest{
#' ds <- cal.dataset(error = "S1", nobs = 50, seed = 1)
#' fit <- cal.bayesian(ds, models = "errors", numSavedSteps = 1000)
#' fit
#' }
#'
#' @import rstan
#'
#' @export

cal.bayesian <- function(calibrationData,
                         models = c("errors", "noerrors", "mixed_errors", "mixed_noerrors"),
                         priors = "Informative",
                         numSavedSteps = 3000,
                         chains = 4,
                         cores = getOption("mc.cores", 1L),
                         MC = NULL) {

  models[models == "mixed"] <- "mixed_noerrors"
  models <- unique(match.arg(models, names(.model_names), several.ok = TRUE))
  .check_cal_data(calibrationData, c("Temperature", "TempError", "D47", "D47error"))
  prior <- .get_priors(priors)

  # The MC argument is kept for older code
  if (!is.null(MC)) {
    warning("MC is deprecated. Use cores instead.")
    if (isTRUE(MC)) cores <- min(chains, parallel::detectCores(), na.rm = TRUE)
  }

  # Mixed models need at least two materials
  mixedModels <- c("mixed_errors", "mixed_noerrors")
  if (any(models %in% mixedModels)) {
    if (length(unique(calibrationData$Material)) < 2) {
      models <- setdiff(models, mixedModels)
      message("Fewer than two materials were found, so the mixed models were not fit.")
    } else if (anyNA(calibrationData$Material)) {
      stop("Material has missing values.")
    }
  }
  if (length(models) == 0) {
    stop("There are no models left to fit.")
  }

  # Data shared by all models. Prior values are ignored when flat = 1.
  flat <- is.null(prior)
  p <- if (flat) list(alpha = c(0, 1), beta = c(0, 1)) else prior
  base <- list(
    N = nrow(calibrationData),
    x = calibrationData$Temperature,
    y = calibrationData$D47,
    flat = as.integer(flat),
    alpha_mu = p$alpha[1],
    alpha_sd = p$alpha[2],
    beta_mu = p$beta[1],
    beta_sd = p$beta[2]
  )
  mats <- levels(factor(calibrationData$Material))

  fits <- list()
  for (m in models) {
    errors <- m %in% c("errors", "mixed_errors")
    mixed <- m %in% mixedModels

    data <- base
    pars <- c("alpha", "beta", "sigma", "log_lik")
    init <- "random"

    if (errors) {
      data$tau <- calibrationData$TempError
      data$y_err <- calibrationData$D47error

      # Start the true temperatures at their measured values
      init <- function() list(x_true = data$x, mu_x = mean(data$x), sigma_x = sd(data$x))
    }
    if (mixed) {
      data$J <- length(mats)
      data$Material <- as.integer(factor(calibrationData$Material))
      pars <- c(pars, "mu_alpha", "mu_beta", "tau_alpha", "tau_beta", "alpha_new", "beta_new")
    }

    fit <- .sample(m, data, pars, numSavedSteps, chains, cores, init = init, adapt_delta = if (mixed) 0.95 else 0.8)

    # Store material labels so reconstructions can match them
    if (mixed) attr(fit, "Material") <- mats
    fits[[.model_names[[m]]]] <- fit
  }

  .finish_fit(fits, prior, chains)
}

#' @rdname cal.bayesian
#' @export

print.bayclump_cal <- function(x, digits = 4, ...) {
  pr <- attr(x, "priors")
  cat("Bayesian clumped isotope calibration\n")
  if (is.list(pr)) {
    cat(sprintf("Priors: alpha ~ N(%g, %g), beta ~ N(%g, %g)\n", pr$alpha[1], pr$alpha[2], pr$beta[1], pr$beta[2]))
  } else {
    cat("Priors:", pr, "\n")
  }

  # Same summary as cal.summary() for each model
  for (m in names(x)) {
    cat("\n", m, "\n", sep = "")
    print(cal.summary(x[[m]]), row.names = FALSE, digits = digits)
  }

  cat("\nConvergence\n")
  print(attr(x, "diagnostics"), row.names = FALSE, digits = digits)
  invisible(x)
}
