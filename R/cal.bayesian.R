#' Bayesian regressions to calibrate the clumped isotopes paleothermometer using
#' \code{stan}.
#'
#' @param calibrationData The target calibration dataset.
#' @param numSavedSteps Number of MCMC iterations to save.
#' @param priors Either \code{Informative}, \code{Weak}, or
#'               \code{Uninformative} on the slope and intercept.
#' @param MC Multicore (TRUE/FALSE)
#'
#' @importFrom loo extract_log_lik relative_eff loo_compare
#' @import parallel
#' @import rstan
#' @importFrom stats sd
#'
#' @export


cal.bayesian <- function(calibrationData,
                         numSavedSteps = 3000,
                         priors = "Informative",
                         MC = TRUE) {

  if(MC){
    options(mc.cores = parallel::detectCores())
  }else{
    if(.Platform$OS.type == "unix"){
      options(mc.cores = 1)
    }
  }

  if (!priors %in% c("Informative", "Weak", "Uninformative")) {
    stop("Priors must be in `Informative`, `Difusse` or `NonInformative`")
  }

  if (priors == "Informative") {
    beta_mu <- 0.039
    beta_sd <- 0.004
    alpha_mu <- 0.231
    alpha_sd <- 0.065
  }

  if (priors == "Uninformative") {
    # Stan actually uses a flat prior
    beta_mu <- 0.01
    beta_sd <- 0.01
    alpha_mu <- 0.01
    alpha_sd <- 0.01
  }

  if (priors == "Weak") {
    beta_mu <- 0.039
    beta_sd <- 1.000
    alpha_mu <- 0.231
    alpha_sd <- 1.000
  }


  ## Models
  fwMod_Errors <- "
  data {
    int<lower=0> N;
    vector[N] y;
    vector[N] x_meas;
    real<lower=0> tau;
    real beta_mu;
    real beta_sd;
    real alpha_mu;
    real alpha_sd;
    real mu_x;
    real sigma_x;
  }

  parameters {
    vector[N] x;
    real alpha;
    real beta;
    real<lower=0> sigma;
  }

  model {
    beta ~ normal(beta_mu, beta_sd);
    alpha ~ normal(alpha_mu, alpha_sd);

    x ~ normal(mu_x, sigma_x);
    x_meas ~ normal(x, tau);
    y ~ normal(alpha + beta * x, sigma);

    sigma ~ cauchy(0, 5);
  }

  generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = normal_lpdf(y[i] | 0, sigma);
    }
  }
"


  fwMod_NE <- "
  data {
    int<lower = 0> N;
    vector[N] x;
    vector[N] y;
    real beta_mu;
    real beta_sd;
    real alpha_mu;
    real alpha_sd;
  }

  parameters {
    real alpha;
    real beta;
    real<lower=0> sigma;
  }

  model {
    alpha ~ normal(alpha_mu, alpha_sd);
    beta ~ normal(beta_mu, beta_sd);
    sigma ~ cauchy(0, 5);
    y ~ normal(alpha + beta * x, sigma);
  }

  generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = normal_lpdf(y[i] | 0, sigma);
    }
  }

"

  fwMod_mixed <- "
  data {
    int<lower=0> N;
    int<lower=0> J;
    vector[N] y;
    vector[N] x;
    array[N] int Material;
    real beta_mu;
    real beta_sd;
    real alpha_mu;
    real alpha_sd;
  }

  parameters {
    real<lower=0> sigma;
    vector[J] alpha;
    vector[J] beta;
  }

  model {
    alpha ~ normal(alpha_mu, alpha_sd);
    beta ~ normal(beta_mu, beta_sd);
    y ~ normal(alpha[Material] + beta[Material].*x, sigma);
    sigma ~ cauchy(0, 5);
  }

  generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = normal_lpdf(y[i] | 0, sigma);
    }
  }
"

  ## Non-Informative

  fwMod_NE2 <- "
  data {
    int<lower = 0> N;
    vector[N] x;
    vector[N] y;
  }

  parameters {
    real alpha;
    real beta;
    real<lower=0> sigma;
  }

  model {
    sigma ~ cauchy(0, 5);
    y ~ normal(alpha + beta * x, sigma);
  }

  generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = normal_lpdf(y[i] | 0, sigma);
    }
  }
"

  fwMod_Errors2 <- "
  data {
    int<lower=0> N;
    vector[N] y;
    vector[N] x_meas;
  }

  parameters {
    vector[N] x;
    real tau;
    real alpha;
    real beta;
    real mu_x;
    real sigma_x;
    real<lower=0> sigma;
  }

  model {
    x ~ normal(mu_x, sigma_x);
    x_meas ~ normal(x, tau);
    y ~ normal(alpha + beta * x, sigma);

    sigma ~ cauchy(0, 5);
  }

  generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = normal_lpdf(y[i] | 0, sigma);
    }
  }
"


  fwMod_mixed2 <- "
  data {
    int<lower=0> N;
    int<lower=0> J;
    vector[N] y;
    vector[N] x;
    array[N] int Material;
  }

  parameters {
    real<lower=0> sigma;
    vector[J] alpha;
    vector[J] beta;
  }

  model {
    y ~ normal(alpha[Material] + beta[Material].*x, sigma);
    sigma ~ cauchy(0, 5);
  }

  generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = normal_lpdf(y[i] | 0, sigma);
    }
  }
"


  # Data

  stan_data_NE <- list(
    N = nrow(calibrationData),
    x = calibrationData$Temperature,
    y = calibrationData$D47,
    beta_mu = beta_mu,
    beta_sd = beta_sd,
    alpha_mu = alpha_mu,
    alpha_sd = alpha_sd
  )

  stan_data_Err <- list(
    N = nrow(calibrationData),
    x_meas = calibrationData$Temperature,
    y = calibrationData$D47,
    tau = mean(calibrationData$TempError),
    beta_mu = beta_mu,
    beta_sd = beta_sd,
    alpha_mu = alpha_mu,
    alpha_sd = alpha_sd,
    mu_x = mean(calibrationData$Temperature),
    sigma_x = sd(calibrationData$Temperature)
  )

  stan_data_mixed <- list(
    N = nrow(calibrationData),
    x = calibrationData$Temperature,
    y = calibrationData$D47,
    J = length(unique(calibrationData$Material)),
    Material = as.numeric(calibrationData$Material),
    beta_mu = beta_mu,
    beta_sd = beta_sd,
    alpha_mu = alpha_mu,
    alpha_sd = alpha_sd
  )

  # Parameters for the run
  nChains <- 2
  burnInSteps <- 1000
  thinSteps <- 1
  nIter <- ceiling(burnInSteps + (numSavedSteps * thinSteps) / nChains)

  # Fit models
  BLM1_E <- stan(
    data = stan_data_Err, model_code = if (priors == "Uninformative") {
      fwMod_Errors2
    } else {
      fwMod_Errors
    },
    chains = nChains, iter = nIter, warmup = burnInSteps,
    thin = thinSteps, pars = c("alpha", "beta", "sigma", "log_lik")
  )

  BLM1_NE <- stan(
    data = stan_data_NE, model_code = if (priors == "Uninformative") {
      fwMod_NE2
    } else {
      fwMod_NE
    },
    chains = nChains, iter = nIter, warmup = burnInSteps,
    thin = thinSteps, pars = c("alpha", "beta", "sigma", "log_lik")
  )

  BLM3 <- stan(
    data = stan_data_mixed, model_code = if (priors == "Uninformative") {
      fwMod_mixed2
    } else {
      fwMod_mixed
    },
    chains = 2, iter = nIter, warmup = burnInSteps,
    thin = thinSteps, pars = c("alpha", "beta", "sigma", "log_lik")
  )

  ## Extract likelihood values for model comparison
  log_lik_BLM1_E <- extract_log_lik(BLM1_E, merge_chains = F)
  r_eff_1_BLM1_E <- relative_eff(log_lik_BLM1_E)
  log_lik_BLM1_NE <- extract_log_lik(BLM1_NE, merge_chains = F)
  r_eff_BLM1_NE <- relative_eff(log_lik_BLM1_NE)
  log_lik_BLM3 <- extract_log_lik(BLM3, merge_chains = F)
  r_eff_BLM3 <- relative_eff(log_lik_BLM3)

  loo_BLM1_E <- loo(log_lik_BLM1_E, r_eff = r_eff_1_BLM1_E)
  loo_BLM1_NE <- loo(log_lik_BLM1_NE, r_eff = r_eff_BLM1_NE)
  loo_BLM3 <- loo(log_lik_BLM3, r_eff = r_eff_BLM3)

  looComp <- loo_compare(list("BLM1_E" = loo_BLM1_E, "BLM1_NE" = loo_BLM1_NE, "BLM3" = loo_BLM3))

  ##return all the relevant objects
  CompleteModelFit <- list(
    "BLM1_fit" = BLM1_E,
    "BLM1_fit_NoErrors" = BLM1_NE,
    "BLM3_fit" = BLM3
  )

  attr(CompleteModelFit, "loo") <- looComp

  return(CompleteModelFit)
}
