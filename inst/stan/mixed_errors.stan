// Mixed model with intercepts and slopes that vary by material, with errors in temperature and D47
data {
  int<lower=0> N;
  int<lower=1> J;                          // Number of materials
  vector[N] x;                             // Measured temperature (10^6/T^2)
  vector<lower=0>[N] tau;                  // 1 SD of x
  vector[N] y;                             // Measured D47
  vector<lower=0>[N] y_err;                // 1 SD of y
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
  vector[N] x_true;
  real mu_x;
  real<lower=0> sigma_x;
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

  // True temperatures and their measurements
  x_true ~ normal(mu_x, sigma_x);
  x ~ normal(x_true, tau);

  // Residual scatter plus the measurement error in D47
  y ~ normal(alpha_c[Material] + beta[Material] .* (x_true - x_mean), sqrt(square(sigma) + square(y_err)));
}

generated quantities {
  vector[N] log_lik;

  // Intercept and slope for a material that is not in the calibration
  real beta_new = normal_rng(mu_beta, tau_beta);
  real alpha_new = normal_rng(mu_alpha_c, tau_alpha) - beta_new * x_mean;

  for (i in 1:N) {
    log_lik[i] = normal_lpdf(y[i] | alpha_c[Material[i]] + beta[Material[i]] * (x_true[i] - x_mean), sqrt(square(sigma) + square(y_err[i])));
  }
}
