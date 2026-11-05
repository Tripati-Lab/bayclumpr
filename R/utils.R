# Internal helpers shared by the calibration and reconstruction functions

# Names of the Bayesian models in the returned objects
.model_names <- c(errors = "BLM1_fit", noerrors = "BLM1_fit_NoErrors", mixed_errors = "BLM3_fit_Errors", mixed_noerrors = "BLM3_fit")

# Models compiled during the session
.models <- new.env(parent = emptyenv())

#' Return a compiled Stan model
#'
#' Uses the models compiled at install time when the package was set up with rstantools. Otherwise it compiles the model from inst/stan once per session.
#'
#' @noRd
.get_model <- function(name) {
  ns <- asNamespace("bayclumpr")
  if (exists("stanmodels", envir = ns, inherits = FALSE)) {
    return(get("stanmodels", envir = ns)[[name]])
  }
  if (is.null(.models[[name]])) {
    file <- system.file("stan", paste0(name, ".stan"), package = "bayclumpr")
    if (!nzchar(file)) {
      stop("The Stan model ", name, ".stan was not found in the installed package. Check that inst/stan in the package source has the six .stan files, then reinstall bayclumpr.")
    }
    # Compiling at run time needs BH, RcppEigen, and StanHeaders, which binary installs of rstan don't include
    .models[[name]] <- tryCatch(
      rstan::stan_model(file = file, model_name = name),
      error = function(e) {
        stop("Compiling the Stan model ", name, " failed: ", conditionMessage(e), " Install bayclumpr from source so the models are compiled at install time, or install the packages BH, RcppEigen, and StanHeaders.", call. = FALSE)
      }
    )
  }
  .models[[name]]
}

#' Return prior means and SDs for alpha and beta, or NULL for flat priors
#'
#' @noRd
.get_priors <- function(priors) {
  if (is.list(priors)) {
    ok <- all(c("alpha", "beta") %in% names(priors)) && length(priors$alpha) == 2 && length(priors$beta) == 2 && priors$alpha[2] > 0 && priors$beta[2] > 0
    if (!isTRUE(ok)) {
      stop("Custom priors must be a list with alpha = c(mean, sd) and beta = c(mean, sd), with positive SDs.")
    }
    return(list(alpha = priors$alpha, beta = priors$beta))
  }
  priors <- match.arg(priors, c("Informative", "Weak", "Flat", "Uninformative"))
  switch(priors,
    Informative = list(alpha = c(0.231, 0.065), beta = c(0.039, 0.004)),
    Weak = list(alpha = c(0.231, 1), beta = c(0.039, 1)),
    NULL
  )
}

#' Check the columns and values of a calibration dataset
#'
#' @noRd
.check_cal_data <- function(data, cols) {
  missing <- setdiff(cols, names(data))
  if (length(missing) > 0) {
    stop("Missing columns in the calibration dataset: ", paste(missing, collapse = ", "), ".")
  }
  if (nrow(data) < 3) {
    stop("The calibration dataset needs at least three observations.")
  }
  if (anyNA(data[, cols])) {
    stop("The calibration dataset has missing values.")
  }
  if (max(data$Temperature) > 20) {
    warning("Temperature should be in 10^6/T^2 (T in K). Values above 20 suggest Celsius or Kelvin. See temp2x().")
  }
  if ("TempError" %in% cols && !all(data$TempError > 0)) {
    stop("TempError must be a positive 1 SD for every observation.")
  }
  if ("D47error" %in% cols && !all(data$D47error >= 0)) {
    stop("D47error must be a non-negative 1 SD for every observation.")
  }
  invisible(TRUE)
}

#' Check the columns and values of a reconstruction dataset
#'
#' @noRd
.check_rec_data <- function(data, cols) {
  missing <- setdiff(cols, names(data))
  if (length(missing) > 0) {
    stop("Missing columns in the reconstruction dataset: ", paste(missing, collapse = ", "), ".")
  }
  if (anyNA(data[, cols])) {
    stop("The reconstruction dataset has missing values.")
  }
  if (!all(data$D47error >= 0)) {
    stop("D47error must be a non-negative 1 SD for every sample.")
  }
  if (!is.null(data$N) && !isTRUE(all(data$N >= 1))) {
    stop("N must be at least 1 for every sample.")
  }
  invisible(TRUE)
}

#' Residual SD of D47 around a calibration line
#'
#' @noRd
.resid_sd <- function(d, alpha, beta) {
  sqrt(sum((d$D47 - alpha - beta * d$Temperature)^2) / (nrow(d) - 2))
}

#' Bootstrap a calibration
#'
#' Resamples the rows of the dataset, fits the model, and returns one row per replicate. Failed replicates are removed with a warning.
#'
#' @noRd
.bootstrap <- function(data, replicates, samples, fit) {
  if (is.null(samples)) samples <- nrow(data)
  reps <- lapply(seq_len(replicates), function(i) {
    d <- data[sample(nrow(data), samples, replace = TRUE), ]
    res <- tryCatch(fit(d), error = function(e) NULL)
    if (is.null(res) || !all(is.finite(res))) NULL else res
  })

  failed <- sum(sapply(reps, is.null))
  if (failed == replicates) {
    stop("All bootstrap replicates failed.")
  }
  if (failed > 0) {
    warning(failed, " of ", replicates, " bootstrap replicates failed and were removed.")
  }
  as.data.frame(do.call(rbind, reps))
}

#' Sample a compiled model with the settings shared by all fits
#'
#' @noRd
.sample <- function(name, data, pars, numSavedSteps, chains, cores, init = "random", adapt_delta = 0.8) {
  warmup <- 1000
  rstan::sampling(
    .get_model(name),
    data = data,
    pars = pars,
    chains = chains,
    cores = cores,
    iter = warmup + ceiling(numSavedSteps / chains),
    warmup = warmup,
    init = init,
    control = list(adapt_delta = adapt_delta),
    refresh = 0
  )
}

#' Summarize convergence for the main parameters of a fit
#'
#' @noRd
.diagnose <- function(fit) {
  pars <- intersect(c("alpha", "beta", "gamma", "sigma"), fit@sim$pars_oi)
  s <- summary(fit, pars = pars)$summary
  data.frame(
    max_Rhat = max(s[, "Rhat"], na.rm = TRUE),
    min_ESS = min(s[, "n_eff"], na.rm = TRUE),
    divergences = rstan::get_num_divergent(fit)
  )
}

#' Add the priors and convergence summary to a list of fits
#'
#' @noRd
.finish_fit <- function(fits, priors, chains) {
  diag <- do.call(rbind, lapply(names(fits), function(m) cbind(model = m, .diagnose(fits[[m]]))))
  bad <- diag$model[diag$max_Rhat > 1.01 | diag$min_ESS < 100 * chains | diag$divergences > 0]
  if (length(bad) > 0) {
    warning("Check convergence for: ", paste(bad, collapse = ", "), ". Print the result for details.")
  }
  attr(fits, "priors") <- if (is.null(priors)) "Flat" else priors
  attr(fits, "diagnostics") <- diag
  class(fits) <- "bayclump_cal"
  fits
}

#' Reconstruct temperatures from calibration draws
#'
#' With a flat prior on 10^6/T^2, the posterior of 10^6/T^2 for each calibration draw is normal, so it is sampled directly. Returns the mean, SD, SE, median, and interval of the temperatures, plus the temperature from the mean slope and intercept.
#'
#' The SE propagates D47error, the standard error of the mean D47. The SD propagates the SD of a single replicate, taken from D47error * sqrt(N) or from D47sd. It is NA without replicate information. The interval is based on the SE.
#'
#' @noRd
.rec_exact <- function(draws, recData, ndraws, level = 0.95, useIon = FALSE, useIonError = FALSE) {
  if (!all(c("alpha", "beta", "sigma") %in% names(draws))) {
    stop("Calibration draws need the columns alpha, beta, and sigma.")
  }
  if (mean(draws$beta) <= 0) {
    stop("The calibration must use Temperature in 10^6/T^2 (T in K).")
  }

  # Mixed models give one set of draws per material
  mixed <- "Material" %in% names(draws)
  if (mixed && (is.null(recData$Material) || !all(as.character(recData$Material) %in% as.character(draws$Material)))) {
    stop("Materials in recData must match those in the calibration. Use newMaterial = TRUE for materials that are not in the calibration.")
  }

  probs <- c(0.5, (1 - level) / 2, 1 - (1 - level) / 2)

  # D47 SD of a single replicate, when replicate information is available
  oneSD <- if (!is.null(recData$N)) {
    recData$D47error * sqrt(recData$N)
  } else if (!is.null(recData$D47sd)) {
    recData$D47sd
  } else {
    rep(NA_real_, nrow(recData))
  }

  res <- t(sapply(seq_len(nrow(recData)), function(i) {
    d <- if (mixed) draws[as.character(draws$Material) == as.character(recData$Material[i]), ] else draws
    k <- sample(nrow(d), ndraws, replace = TRUE)
    mu <- recData$D47[i] - d$alpha[k]

    # Remove the effect of the covariate, sampling its error when needed
    if (useIon) {
      ion <- if (useIonError) rnorm(ndraws, recData$Ion[i], recData$IonError[i]) else recData$Ion[i]
      mu <- mu - d$gamma[k] * ion
    }

    # Temperatures for a given measurement error in D47, keeping physically possible values of 10^6/T^2
    temps <- function(err) {
      x <- (mu + rnorm(ndraws, 0, sqrt(d$sigma[k]^2 + err^2))) / d$beta[k]
      sqrt(10^6 / x[x > 0]) - 273.15
    }
    tempSE <- temps(recData$D47error[i])
    tempSD <- if (is.na(oneSD[i])) NA else sd(temps(oneSD[i]))

    # Temperature from the mean slope and intercept
    plug <- sqrt(10^6 * mean(d$beta) / (recData$D47[i] - mean(d$alpha))) - 273.15

    c(mean(tempSE), tempSD, sd(tempSE), quantile(tempSE, probs, names = FALSE), plug, ndraws - length(tempSE))
  }))

  if (sum(res[, 8]) > 0) {
    warning(sum(res[, 8]), " draws had 10^6/T^2 at or below zero and were removed.")
  }

  out <- cbind.data.frame(
    Sample = recData$Sample,
    D47 = recData$D47,
    D47error = recData$D47error,
    meanTemp = res[, 1],
    SD = res[, 2],
    SE = res[, 3],
    median = res[, 4],
    lower = res[, 5],
    upper = res[, 6]
  )
  attr(out, "plug") <- res[, 7]
  out
}
