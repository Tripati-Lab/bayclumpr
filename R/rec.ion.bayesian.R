#' This function generates temperature predictions (in 10^6/T2) based on a
#' calibration dataset and target D47. It includes measured error in the target D47
#' and ion as an additional predictor.
#'
#' @param calModel The fitted Stan model from calibration.
#' @param recData The reconstruction dataset.
#' @param iter Number of sampling iterations.
#' @param postcalsamples Number of posterior calibration samples to use.
#' @param MC Enable multicore (default TRUE).
#'
#' @import rstan
#' @import parallel
#'
#' @export
rec.ion.bayesian <- function(calModel,
                         recData,
                         iter = 1000,
                         postcalsamples = NULL,
                         MC = TRUE) {

  if (MC) {
    options(mc.cores = parallel::detectCores())
  } else {
    options(mc.cores = 1)
  }

  vects.params <- extract(calModel)

  Model <- "
  data {
    int<lower=0> n;
    vector[n] y_mes;
    vector[n] y_err;
    vector[n] ion_meas;
    vector[n] ion_err;
    int<lower=0> posts;
    vector[posts] alpha;
    vector[posts] beta;
    vector[posts] gamma;
    vector[posts] sigma;
  }

  parameters {
    matrix[n, posts] x_new;
    matrix[n, posts] ion_new;
  }

  model {
    vector[posts] y_hat;
    for (i in 1:n) {
      // Observation model includes prediction and error
      y_hat = alpha + beta .* x_new[i,]' + gamma .* ion_new[i,]';
      y_mes[i] ~ normal(y_hat, sqrt(sigma^2 + y_err[i]^2));
      ion_meas[i] ~ normal(ion_new[i,]', ion_err[i]);
    }
  }
  "

  totsamp <- length(vects.params[[1]])
  seqSamples <- if (!is.null(postcalsamples)) {
    sample(1:totsamp, postcalsamples, replace = FALSE)
  } else {
    1:totsamp
  }


    stan_data <- list(
      n = nrow(recData),
      y_mes = recData$D47,
      y_err = recData$D47error,
      ion_meas = recData$Ion,
      ion_err = recData$IonError,
      posts = length(seqSamples),
      alpha = vects.params$alpha[seqSamples],
      beta = vects.params$beta[seqSamples],
      gamma = vects.params$gamma[seqSamples],
      sigma = vects.params$sigma[seqSamples]
    )

    fit <- stan(
      data = stan_data,
      model_code = Model,
      chains = 2,
      iter = iter,
      warmup = floor(iter / 2),
      control = list(adapt_delta = 0.9)
    )

    samples <- extract(fit)
    Xouts <- samples$x_new

    temp_mat <- sqrt(1e6 / Xouts) - 273.15

    temp_summary <- t(apply(temp_mat, 2, function(x) c(mean = mean(x), sd = sd(x))))

    return(data.frame(
      Sample = recData$Sample,
      D47 = recData$D47,
      D47error = recData$D47error,
      meanTemp = temp_summary[, "mean"],
      error = temp_summary[, "sd"]
    ))
}
