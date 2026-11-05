// Linear calibration with errors in temperature and D47
data {
  int<lower=0> N;
  vector[N] x;                // Measured temperature (10^6/T^2)
  vector<lower=0>[N] tau;     // 1 SD of x
  vector[N] y;                // Measured D47
  vector<lower=0>[N] y_err;   // 1 SD of y
  int<lower=0, upper=1> flat; // 1 for flat priors on alpha and beta
  real alpha_mu;
  real<lower=0> alpha_sd;
  real beta_mu;
  real<lower=0> beta_sd;
}

parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
  vector[N] x_true;
  real mu_x;
  real<lower=0> sigma_x;
}

model {
  if (flat == 0) {
    alpha ~ normal(alpha_mu, alpha_sd);
    beta ~ normal(beta_mu, beta_sd);
  }
  sigma ~ cauchy(0, 5);

  // True temperatures and their measurements
  x_true ~ normal(mu_x, sigma_x);
  x ~ normal(x_true, tau);

  // Residual scatter plus the measurement error in D47
  y ~ normal(alpha + beta * x_true, sqrt(square(sigma) + square(y_err)));
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = normal_lpdf(y[i] | alpha + beta * x_true[i], sqrt(square(sigma) + square(y_err[i])));
  }
}
