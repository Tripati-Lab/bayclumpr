#' Bayesian regressions to calibrate the clumped isotopes paleothermometer using
#' \code{stan}.
#'
#' @param calibrationData The target calibration dataset. It needs the columns \code{Temperature} (10^6/T^2, with T in K), \code{TempError} (1 SD of \code{Temperature}), \code{D47}, \code{D47error} (1 SD of \code{D47}), and \code{Ion}.
#' @param numSavedSteps Number of MCMC iterations to save.
#' @param MC Multicore (TRUE/FALSE)
#' @param IonError Ion measurement error (1 SD). Use a single value for all samples or one value per sample. Use \code{NULL} to treat \code{Ion} as known without error.
#'
#' @details Both models account for errors in \code{Temperature} and \code{D47}. Their \code{sigma} is the scatter left after removing the measurement error in \code{D47}.
#'
#' @import parallel
#' @import rstan
#' @importFrom stats sd
#'
#' @export


cal.ion.bayesian <- function(calibrationData,
                             numSavedSteps = 3000,
                             MC = TRUE,
                             IonError = NULL) {
  if(MC){
    options(mc.cores = parallel::detectCores())
  } else {
    if(.Platform$OS.type == "unix"){
      options(mc.cores = 1)
    }
  }

  # Temperature errors are used per observation and must be positive
  if (!isTRUE(all(calibrationData$TempError > 0))) {
    stop("TempError must be a positive 1 SD for every observation.")
  }
  if (!isTRUE(all(calibrationData$D47error >= 0))) {
    stop("D47error must be a non-negative 1 SD for every observation.")
  }

  # Flat prior
  beta_mu <- 0.01
  beta_sd <- 0.01
  alpha_mu <- 0.01
  alpha_sd <- 0.01

  # Model with ion error
  fwMod_IonErr <- "
  data {
    int<lower=0> N;
    vector[N] y;
    vector<lower=0>[N] y_err;
    vector[N] x_meas;
    vector[N] ion_meas;
    vector<lower=0>[N] tau;
    vector<lower=0>[N] ion_tau;
    real mu_x;
    real<lower=0> sigma_x;
    real mu_ion;
    real<lower=0> sigma_ion;
  }
  parameters {
    vector[N] x;
    vector[N] ion;
    real alpha;
    real beta;
    real gamma;
    real<lower=0> sigma;
  }
  model {
    x ~ normal(mu_x, sigma_x);
    ion ~ normal(mu_ion, sigma_ion);
    x_meas ~ normal(x, tau);
    ion_meas ~ normal(ion, ion_tau);
    y ~ normal(alpha + beta * x + gamma * ion, sqrt(square(sigma) + square(y_err)));
    sigma ~ cauchy(0, 5);
  }
  generated quantities {
    vector[N] log_lik;
    for (i in 1:N) {
      log_lik[i] = normal_lpdf(y[i] | alpha + beta * x[i] + gamma * ion[i], sqrt(square(sigma) + square(y_err[i])));
    }
  }
  "

  # Model without ion error
  fwMod_NoIonErr <- "
  data {
    int<lower=0> N;
    vector[N] y;
    vector<lower=0>[N] y_err;
    vector[N] x_meas;
    vector[N] ion;
    vector<lower=0>[N] tau;
    real mu_x;
    real<lower=0> sigma_x;
  }
  parameters {
    vector[N] x;
    real alpha;
    real beta;
    real gamma;
    real<lower=0> sigma;
  }
  model {
    x ~ normal(mu_x, sigma_x);
    x_meas ~ normal(x, tau);
    y ~ normal(alpha + beta * x + gamma * ion, sqrt(square(sigma) + square(y_err)));
    sigma ~ cauchy(0, 5);
  }
  generated quantities {
    vector[N] log_lik;
    for (i in 1:N) {
      log_lik[i] = normal_lpdf(y[i] | alpha + beta * x[i] + gamma * ion[i], sqrt(square(sigma) + square(y_err[i])));
    }
  }
  "

  # Determine which model to use based on whether IonError is provided
  if (!is.null(IonError)) {
    # Use one error per sample, or the same error for all samples
    if (!length(IonError) %in% c(1, nrow(calibrationData)) || !all(IonError > 0)) {
      stop("IonError must be positive, with length 1 or one value per observation.")
    }
    stan_data <- list(
      N = nrow(calibrationData),
      x_meas = calibrationData$Temperature,
      ion_meas = calibrationData$Ion,
      y = calibrationData$D47,
      y_err = calibrationData$D47error,
      tau = calibrationData$TempError,
      ion_tau = rep(IonError, length.out = nrow(calibrationData)),
      mu_x = mean(calibrationData$Temperature),
      sigma_x = sd(calibrationData$Temperature),
      mu_ion = mean(calibrationData$Ion),
      sigma_ion = sd(calibrationData$Ion)
    )
    model_code <- fwMod_IonErr
  } else {
    # No ion error
    stan_data <- list(
      N = nrow(calibrationData),
      x_meas = calibrationData$Temperature,
      ion = calibrationData$Ion,
      y = calibrationData$D47,
      y_err = calibrationData$D47error,
      tau = calibrationData$TempError,
      mu_x = mean(calibrationData$Temperature),
      sigma_x = sd(calibrationData$Temperature)
    )
    model_code <- fwMod_NoIonErr
  }

  # Parameters for the run
  nChains <- 2
  burnInSteps <- 1000
  thinSteps <- 1
  nIter <- ceiling(burnInSteps + (numSavedSteps * thinSteps) / nChains)

  # Fit models
  BLM1_fit <- stan(
    data = stan_data, model_code = model_code,
    chains = nChains, iter = nIter, warmup = burnInSteps,
    thin = thinSteps, pars = c("alpha", "beta", "gamma", "sigma")
  )

  # Return all the relevant objects
  CompleteModelFit <- list(
    "BLM1_fit" = BLM1_fit
  )
  return(CompleteModelFit)
}
