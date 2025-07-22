#' Bayesian regressions to calibrate the clumped isotopes paleothermometer using
#' \code{stan}.
#'
#' @param calibrationData The target calibration dataset.
#' @param numSavedSteps Number of MCMC iterations to save.
#' @param MC Multicore (TRUE/FALSE)
#' @param useIonError Whether or not to account for ion measurement error.
#'
#' @import parallel
#' @import rstan
#' @importFrom stats sd
#'
#' @export


cal.ion.bayesian <- function(calibrationData,
                             numSavedSteps = 3000,
                             MC = TRUE,
                             useIonError = FALSE) {

  if(MC){
    options(mc.cores = parallel::detectCores())
  } else {
    if(.Platform$OS.type == "unix"){
      options(mc.cores = 1)
    }
  }

  # Flat prior
  beta_mu <- 0.01
  beta_sd <- 0.01
  alpha_mu <- 0.01
  alpha_sd <- 0.01

  fwMod_IonErr <- "
  data {
    int<lower=0> N;
    vector[N] y;
    vector[N] x_meas;
    vector[N] ion_meas;
    real<lower=0> tau;
    real<lower=0> ion_tau;
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
    y ~ normal(alpha + beta * x + gamma * ion, sigma);
    sigma ~ cauchy(0, 5);
  }

  generated quantities {
    vector[N] log_lik;
    for (i in 1:N) {
      log_lik[i] = normal_lpdf(y[i] | alpha + beta * x[i] + gamma * ion[i], sigma);
    }
  }
  "

  fwMod_NoIonErr <- "
  data {
    int<lower=0> N;
    vector[N] y;
    vector[N] x_meas;
    vector[N] ion;
    real<lower=0> tau;
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
    y ~ normal(alpha + beta * x + gamma * ion, sigma);
    sigma ~ cauchy(0, 5);
  }

  generated quantities {
    vector[N] log_lik;
    for (i in 1:N) {
      log_lik[i] = normal_lpdf(y[i] | alpha + beta * x[i] + gamma * ion[i], sigma);
    }
  }
  "

  if (useIonError) {
    stan_data <- list(
      N = nrow(calibrationData),
      x_meas = calibrationData$Temperature,
      ion_meas = calibrationData$Ion,
      y = calibrationData$D47,
      tau = sd(calibrationData$TempError),
      ion_tau = sd(calibrationData$IonError),
      mu_x = mean(calibrationData$Temperature),
      sigma_x = sd(calibrationData$Temperature),
      mu_ion = mean(calibrationData$Ion),
      sigma_ion = sd(calibrationData$Ion)
    )
    model_code <- fwMod_IonErr
  } else {
    stan_data <- list(
      N = nrow(calibrationData),
      x_meas = calibrationData$Temperature,
      ion = calibrationData$Ion,
      y = calibrationData$D47,
      tau = sd(calibrationData$TempError),
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

  ## Return all the relevant objects
  CompleteModelFit <- list(
    "BLM1_fit" = BLM1_fit
  )

  return(CompleteModelFit)
}
