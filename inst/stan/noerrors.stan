// Linear calibration that ignores measurement errors
data {
  int<lower=0> N;
  vector[N] x;                // Temperature (10^6/T^2)
  vector[N] y;                // D47
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
  y ~ normal(alpha_c + beta * (x - x_mean), sigma);
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = normal_lpdf(y[i] | alpha_c + beta * (x[i] - x_mean), sigma);
  }
}
