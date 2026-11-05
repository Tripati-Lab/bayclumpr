// Mixed model with intercepts and slopes that vary by material
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

parameters {
  real mu_alpha;
  real mu_beta;
  real<lower=0> tau_alpha;
  real<lower=0> tau_beta;
  vector[J] z_alpha;
  vector[J] z_beta;
  real<lower=0> sigma;
}

transformed parameters {
  // Non-centered parameterization
  vector[J] alpha = mu_alpha + tau_alpha * z_alpha;
  vector[J] beta = mu_beta + tau_beta * z_beta;
}

model {
  if (flat == 0) {
    mu_alpha ~ normal(alpha_mu, alpha_sd);
    mu_beta ~ normal(beta_mu, beta_sd);
  }
  tau_alpha ~ normal(0, 0.1);
  tau_beta ~ normal(0, 0.01);
  z_alpha ~ std_normal();
  z_beta ~ std_normal();
  sigma ~ cauchy(0, 5);
  y ~ normal(alpha[Material] + beta[Material] .* x, sigma);
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = normal_lpdf(y[i] | alpha[Material[i]] + beta[Material[i]] * x[i], sigma);
  }
}
