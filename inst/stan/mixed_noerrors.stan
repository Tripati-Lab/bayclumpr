// Mixed model with intercepts and slopes that vary by material, ignoring measurement errors
data {
  int<lower=0> N;
  int<lower=1> J;                          // Number of materials
  vector[N] x;                             // Temperature (10^6/T^2)
  vector[N] y;                             // D47
  array[N] int<lower=1, upper=J> Material; // Material of each observation
  int<lower=0, upper=1> flat;              // 1 for flat priors on the shared means
  real alpha_mu;
  real<lower=0> alpha_sd;
  real beta_mu;
  real<lower=0> beta_sd;
}

transformed data {
  // Temperature is centered to improve sampling
  real x_mean = mean(x);
}

parameters {
  real mu_alpha_c;
  real mu_beta;
  real<lower=0> tau_alpha;
  real<lower=0> tau_beta;
  vector[J] z_alpha;
  vector[J] z_beta;
  real<lower=0> sigma;
}

transformed parameters {
  // Non-centered parameterization. Intercepts vary at the mean temperature.
  vector[J] alpha_c = mu_alpha_c + tau_alpha * z_alpha;
  vector[J] beta = mu_beta + tau_beta * z_beta;

  // Intercepts at 10^6/T^2 = 0, where the priors are set
  vector[J] alpha = alpha_c - beta * x_mean;
  real mu_alpha = mu_alpha_c - mu_beta * x_mean;
}

model {
  if (flat == 0) {
    target += normal_lpdf(mu_alpha | alpha_mu, alpha_sd);
    target += normal_lpdf(mu_beta | beta_mu, beta_sd);
  }
  tau_alpha ~ normal(0, 0.1);
  tau_beta ~ normal(0, 0.01);
  z_alpha ~ std_normal();
  z_beta ~ std_normal();
  sigma ~ cauchy(0, 5);
  y ~ normal(alpha_c[Material] + beta[Material] .* (x - x_mean), sigma);
}

generated quantities {
  vector[N] log_lik;

  // Intercept and slope for a material that is not in the calibration
  real beta_new = normal_rng(mu_beta, tau_beta);
  real alpha_new = normal_rng(mu_alpha_c, tau_alpha) - beta_new * x_mean;

  for (i in 1:N) {
    log_lik[i] = normal_lpdf(y[i] | alpha_c[Material[i]] + beta[Material[i]] * (x[i] - x_mean), sigma);
  }
}
