#' Bayesian calibration of the clumped isotope paleothermometer
#'
#' Fits up to three Bayesian linear regressions of D47 on temperature using Stan.
#'
#' @param calibrationData The calibration dataset. It needs the columns \code{Temperature} (10^6/T^2, with T in K), \code{TempError} (1 SD of \code{Temperature}), \code{D47}, and \code{D47error} (1 SD of \code{D47}). The mixed model also needs \code{Material}.
#' @param models The models to fit. Use any of \code{"errors"}, \code{"noerrors"}, and \code{"mixed"}.
#' @param priors Priors on the intercept (alpha) and slope (beta). Use \code{"Informative"}, \code{"Weak"}, or \code{"Flat"}, or a list such as \code{list(alpha = c(0.231, 0.065), beta = c(0.039, 0.004))} with the mean and SD of each normal prior. \code{"Uninformative"} is kept as another name for \code{"Flat"}.
#' @param numSavedSteps Total number of posterior draws to keep across all chains.
#' @param chains Number of MCMC chains.
#' @param cores Number of cores used to run the chains in parallel.
#' @param MC Deprecated. Use \code{cores} instead.
#' @param x A \code{bayclump_cal} object.
#' @param digits Number of decimal places to print.
#' @param ... Not used.
#'
#' @details The three models are:
#' \itemize{
#'   \item \code{"errors"} (\code{BLM1_fit}) accounts for errors in both \code{Temperature} and \code{D47}. Its \code{sigma} is the scatter left after removing the measurement error in \code{D47}. The true temperatures have a normal prior with a mean and SD estimated from the data.
#'   \item \code{"noerrors"} (\code{BLM1_fit_NoErrors}) ignores both errors, so its \code{sigma} includes the measurement error.
#'   \item \code{"mixed"} (\code{BLM3_fit}) lets intercepts and slopes vary by \code{Material} around shared means (partial pooling). It ignores measurement errors. The priors apply to the shared means. The SDs among materials have half-normal priors with scales of 0.1 (intercepts) and 0.01 (slopes). With fewer than five materials, these SDs are mostly set by their priors. The model is skipped when there are fewer than two materials.
#' }
#'
#' The \code{"Informative"} values are the ones used in earlier versions of this package. Check that they suit the reference frame of your data, or set your own.
#'
#' A warning is given when a model shows signs of poor convergence: Rhat above 1.01, fewer than 100 effective draws per chain, or divergent transitions.
#'
#' @return An object of class \code{bayclump_cal}. It is a list of \code{stanfit} objects named \code{BLM1_fit}, \code{BLM1_fit_NoErrors}, and \code{BLM3_fit}, with the priors and a convergence summary stored as attributes. Use \code{cal.draws()} to get the posterior draws and \code{cal.compare()} to compare models.
#'
#' @import rstan
#'
#' @export

cal.bayesian <- function(calibrationData,
                         models = c("errors", "noerrors", "mixed"),
                         priors = "Informative",
                         numSavedSteps = 3000,
                         chains = 4,
                         cores = getOption("mc.cores", 1L),
                         MC = NULL) {

  models <- match.arg(models, several.ok = TRUE)
  .check_cal_data(calibrationData, c("Temperature", "TempError", "D47", "D47error"))
  prior <- .get_priors(priors)

  # The MC argument is kept for older code
  if (!is.null(MC)) {
    warning("MC is deprecated. Use cores instead.")
    if (isTRUE(MC)) cores <- min(chains, parallel::detectCores(), na.rm = TRUE)
  }

  # A mixed model needs at least two materials
  if ("mixed" %in% models && length(unique(calibrationData$Material)) < 2) {
    models <- setdiff(models, "mixed")
    message("Fewer than two materials were found, so the mixed model was not fit.")
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

  fits <- list()

  if ("errors" %in% models) {
    data <- c(base, list(tau = calibrationData$TempError, y_err = calibrationData$D47error))

    # Start the true temperatures at their measured values
    init <- function() list(x_true = data$x, mu_x = mean(data$x), sigma_x = sd(data$x))
    fits$BLM1_fit <- .sample("errors", data, c("alpha", "beta", "sigma", "log_lik"), numSavedSteps, chains, cores, init = init)
  }

  if ("noerrors" %in% models) {
    fits$BLM1_fit_NoErrors <- .sample("noerrors", base, c("alpha", "beta", "sigma", "log_lik"), numSavedSteps, chains, cores)
  }

  if ("mixed" %in% models) {
    if (anyNA(calibrationData$Material)) {
      stop("Material has missing values.")
    }
    mats <- levels(factor(calibrationData$Material))
    data <- c(base, list(J = length(mats), Material = as.integer(factor(calibrationData$Material))))
    pars <- c("alpha", "beta", "mu_alpha", "mu_beta", "tau_alpha", "tau_beta", "sigma", "log_lik")
    fit <- .sample("mixed", data, pars, numSavedSteps, chains, cores, adapt_delta = 0.95)

    # Store material labels so reconstructions can match them
    attr(fit, "Material") <- mats
    fits$BLM3_fit <- fit
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

  # Posterior summary of the main parameters in each model
  for (m in names(x)) {
    pars <- intersect(c("alpha", "beta", "gamma", "sigma"), x[[m]]@sim$pars_oi)
    s <- summary(x[[m]], pars = pars, probs = c(0.025, 0.975))$summary
    mats <- attr(x[[m]], "Material")
    if (!is.null(mats)) {
      cat("\n", m, " (materials: ", paste(mats, collapse = ", "), ")\n", sep = "")
    } else {
      cat("\n", m, "\n", sep = "")
    }
    print(round(s[, c("mean", "sd", "2.5%", "97.5%")], digits))
  }

  cat("\nConvergence\n")
  print(attr(x, "diagnostics"), row.names = FALSE, digits = digits)
  invisible(x)
}
