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

transformed data {
  // Temperature is centered to improve sampling
  real x_mean = mean(x);
}

parameters {
  real alpha_c;
  real beta;
  real<lower=0> sigma;
  vector[N] x_true;
  real mu_x;
  real<lower=0> sigma_x;
}

transformed parameters {
  // Intercept at 10^6/T^2 = 0, where the priors are set
  real alpha = alpha_c - beta * x_mean;
}

model {
  if (flat == 0) {
    target += normal_lpdf(alpha | alpha_mu, alpha_sd);
    target += normal_lpdf(beta | beta_mu, beta_sd);
  }
  sigma ~ cauchy(0, 5);

  // True temperatures and their measurements
  x_true ~ normal(mu_x, sigma_x);
  x ~ normal(x_true, tau);

  // Residual scatter plus the measurement error in D47
  y ~ normal(alpha_c + beta * (x_true - x_mean), sqrt(square(sigma) + square(y_err)));
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = normal_lpdf(y[i] | alpha_c + beta * (x_true[i] - x_mean), sqrt(square(sigma) + square(y_err[i])));
  }
}
