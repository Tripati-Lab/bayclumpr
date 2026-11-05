#' Bayesian regressions to calibrate the clumped isotopes paleothermometer using
#' \code{stan}.
#'
#' @param calibrationData The target calibration dataset. It needs the columns \code{Temperature} (10^6/T^2, with T in K), \code{TempError} (1 SD of \code{Temperature}), \code{D47}, \code{D47error} (1 SD of \code{D47}), and \code{Material}.
#' @param numSavedSteps Number of MCMC iterations to save.
#' @param priors Either \code{Informative}, \code{Weak}, or
#'               \code{Uninformative} on the slope and intercept.
#' @param MC Multicore (TRUE/FALSE)
#'
#' @details \code{BLM1_fit} accounts for errors in both \code{Temperature} and \code{D47}. Its \code{sigma} is the scatter left after removing the measurement error in \code{D47}. \code{BLM1_fit_NoErrors} ignores both errors, so its \code{sigma} includes the measurement error. \code{BLM3_fit} is a mixed model with intercepts and slopes that vary by \code{Material} around a shared mean (partial pooling). The priors set by \code{priors} apply to these shared means. The SDs among materials have half-normal priors with scales of 0.1 (intercepts) and 0.01 (slopes). With fewer than five materials, these SDs are mostly set by their priors. \code{BLM3_fit} is only fit when there are at least two materials. Otherwise it is \code{NULL}.
#'
#' The \code{loo} attribute compares the fitted models using PSIS-LOO. In \code{BLM1_fit}, the pointwise log-likelihood is conditional on the latent true temperatures. A warning is given when any Pareto k value is above 0.7, since the comparison may then be unreliable.
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
    stop("Priors must be in `Informative`, `Weak`, or `Uninformative`")
  }

  # Temperature errors are used per observation and must be positive
  if (!isTRUE(all(calibrationData$TempError > 0))) {
    stop("TempError must be a positive 1 SD for every observation.")
  }
  if (!isTRUE(all(calibrationData$D47error >= 0))) {
    stop("D47error must be a non-negative 1 SD for every observation.")
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


  # Models
  fwMod_Errors <- "
  data {
    int<lower=0> N;
    vector[N] y;
    vector<lower=0>[N] y_err;
    vector[N] x_meas;
    vector<lower=0>[N] tau;
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
    y ~ normal(alpha + beta * x, sqrt(square(sigma) + square(y_err)));

    sigma ~ cauchy(0, 5);
  }

  generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = normal_lpdf(y[i] | alpha + beta * x[i], sqrt(square(sigma) + square(y_err[i])));
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
    log_lik[i] = normal_lpdf(y[i] | alpha + beta * x[i], sigma);
    }
  }

"

  fwMod_mixed <- "
  data {
    int<lower=0> N;
    int<lower=1> J;
    vector[N] y;
    vector[N] x;
    array[N] int<lower=1, upper=J> Material;
    real beta_mu;
    real beta_sd;
    real alpha_mu;
    real alpha_sd;
  }

  parameters {
    real<lower=0> sigma;
    real mu_alpha;
    real mu_beta;
    real<lower=0> tau_alpha;
    real<lower=0> tau_beta;
    vector[J] z_alpha;
    vector[J] z_beta;
  }

  transformed parameters {
    vector[J] alpha = mu_alpha + tau_alpha * z_alpha;
    vector[J] beta = mu_beta + tau_beta * z_beta;
  }

  model {
    mu_alpha ~ normal(alpha_mu, alpha_sd);
    mu_beta ~ normal(beta_mu, beta_sd);
    tau_alpha ~ normal(0, 0.1);
    tau_beta ~ normal(0, 0.01);
    z_alpha ~ std_normal();
    z_beta ~ std_normal();
    y ~ normal(alpha[Material] + beta[Material].*x, sigma);
    sigma ~ cauchy(0, 5);
  }

  generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = normal_lpdf(y[i] | alpha[Material[i]] + beta[Material[i]] * x[i], sigma);
    }
  }
"

  # Non-Informative

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
    log_lik[i] = normal_lpdf(y[i] | alpha + beta * x[i], sigma);
    }
  }
"

  fwMod_Errors2 <- "
  data {
    int<lower=0> N;
    vector[N] y;
    vector<lower=0>[N] y_err;
    vector[N] x_meas;
    vector<lower=0>[N] tau;
  }

  parameters {
    vector[N] x;
    real alpha;
    real beta;
    real mu_x;
    real<lower=0> sigma_x;
    real<lower=0> sigma;
  }

  model {
    x ~ normal(mu_x, sigma_x);
    x_meas ~ normal(x, tau);
    y ~ normal(alpha + beta * x, sqrt(square(sigma) + square(y_err)));

    sigma ~ cauchy(0, 5);
  }

  generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = normal_lpdf(y[i] | alpha + beta * x[i], sqrt(square(sigma) + square(y_err[i])));
    }
  }
"


  fwMod_mixed2 <- "
  data {
    int<lower=0> N;
    int<lower=1> J;
    vector[N] y;
    vector[N] x;
    array[N] int<lower=1, upper=J> Material;
  }

  parameters {
    real<lower=0> sigma;
    real mu_alpha;
    real mu_beta;
    real<lower=0> tau_alpha;
    real<lower=0> tau_beta;
    vector[J] z_alpha;
    vector[J] z_beta;
  }

  transformed parameters {
    vector[J] alpha = mu_alpha + tau_alpha * z_alpha;
    vector[J] beta = mu_beta + tau_beta * z_beta;
  }

  model {
    tau_alpha ~ normal(0, 0.1);
    tau_beta ~ normal(0, 0.01);
    z_alpha ~ std_normal();
    z_beta ~ std_normal();
    y ~ normal(alpha[Material] + beta[Material].*x, sigma);
    sigma ~ cauchy(0, 5);
  }

  generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = normal_lpdf(y[i] | alpha[Material[i]] + beta[Material[i]] * x[i], sigma);
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
    y_err = calibrationData$D47error,
    tau = calibrationData$TempError,
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
    Material = as.integer(factor(calibrationData$Material)),
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

  # A mixed model needs at least two materials
  if (stan_data_mixed$J > 1) {
    BLM3 <- stan(
      data = stan_data_mixed, model_code = if (priors == "Uninformative") {
        fwMod_mixed2
      } else {
        fwMod_mixed
      },
      chains = 2, iter = nIter, warmup = burnInSteps,
      thin = thinSteps, control = list(adapt_delta = 0.95),
      pars = c("alpha", "beta", "mu_alpha", "mu_beta", "tau_alpha", "tau_beta", "sigma", "log_lik")
    )

    # Store material labels so rec.bayesian() can match them
    attr(BLM3, "Material") <- levels(factor(calibrationData$Material))
  } else {
    BLM3 <- NULL
    message("Fewer than two materials were found, so the mixed model (BLM3) was not fit.")
  }

  # Compare the fitted models using PSIS-LOO
  fits <- list("BLM1_E" = BLM1_E, "BLM1_NE" = BLM1_NE, "BLM3" = BLM3)
  fits <- fits[!sapply(fits, is.null)]
  loos <- lapply(fits, function(fit) {
    log_lik <- extract_log_lik(fit, merge_chains = FALSE)
    loo(log_lik, r_eff = relative_eff(log_lik))
  })
  looComp <- loo_compare(loos)

  # Flag unreliable PSIS-LOO estimates
  k <- unlist(lapply(loos, loo::pareto_k_values))
  if (any(k > 0.7)) {
    warning("Some Pareto k values are above 0.7, so the LOO comparison may be unreliable.")
  }

  # Return all the relevant objects
  CompleteModelFit <- list(
    "BLM1_fit" = BLM1_E,
    "BLM1_fit_NoErrors" = BLM1_NE,
    "BLM3_fit" = BLM3
  )

  attr(CompleteModelFit, "loo") <- looComp

  return(CompleteModelFit)
}
