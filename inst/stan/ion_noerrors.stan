// Calibration with a second predictor (Ion), with errors in temperature and D47 but not in Ion
data {
  int<lower=0> N;
  vector[N] x;                  // Measured temperature (10^6/T^2)
  vector<lower=0>[N] tau;       // 1 SD of x
  vector[N] ion;                // Ion, treated as known
  vector[N] y;                  // Measured D47
  vector<lower=0>[N] y_err;     // 1 SD of y
}

transformed data {
  // Both predictors are centered to improve sampling
  real x_mean = mean(x);
  real ion_mean = mean(ion);
}

parameters {
  real alpha_c;
  real beta;
  real gamma;
  real<lower=0> sigma;
  vector[N] x_true;
  real mu_x;
  real<lower=0> sigma_x;
}

transformed parameters {
  // Intercept at 10^6/T^2 = 0 and Ion = 0
  real alpha = alpha_c - beta * x_mean - gamma * ion_mean;
}

model {
  sigma ~ cauchy(0, 5);

  // True temperatures and their measurements
  x_true ~ normal(mu_x, sigma_x);
  x ~ normal(x_true, tau);

  // Residual scatter plus the measurement error in D47
  y ~ normal(alpha_c + beta * (x_true - x_mean) + gamma * (ion - ion_mean), sqrt(square(sigma) + square(y_err)));
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = normal_lpdf(y[i] | alpha_c + beta * (x_true[i] - x_mean) + gamma * (ion[i] - ion_mean), sqrt(square(sigma) + square(y_err[i])));
  }
}
