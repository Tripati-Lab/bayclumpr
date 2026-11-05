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

parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}

model {
  if (flat == 0) {
    alpha ~ normal(alpha_mu, alpha_sd);
    beta ~ normal(beta_mu, beta_sd);
  }
  sigma ~ cauchy(0, 5);
  y ~ normal(alpha + beta * x, sigma);
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = normal_lpdf(y[i] | alpha + beta * x[i], sigma);
  }
}
