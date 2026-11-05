#' This function generates temperature predictions (in degrees C) based on a
#' calibration dataset and target D47. Note that this approach additionally
#' accounts for measured error in the target D47. This approach is congruent
#' with the one used in McClelland et al. (2022).
#'
#' @param calModel The stan model to be analyzed. The calibration must use \code{Temperature} in 10^6/T^2 (T in K).
#' @param recData The reconstruction dataset. It needs the columns \code{Sample}, \code{D47}, and \code{D47error} (1 SD of the D47 value; if D47 is a mean of replicates, use its standard error). Mixed models also need \code{Material}, with the same labels used in the calibration.
#' @param iter Number of replicates to retain.
#' @param mixed whether the model \code{calModel} is mixed or not.
#' @param postcalsamples Number of posterior samples to analyze from the calibration step.
#' @param MC Multicore (TRUE/FALSE)
#'
#' @return A \code{data.frame} with \code{meanTemp} (posterior mean, in degrees C) and \code{error} (1 SD) for each sample. The error includes uncertainty in the calibration parameters, the residual SD of the calibration (\code{sigma}), and \code{D47error}.
#'
#' @import rstan
#' @import parallel
#'
#' @export

rec.bayesian <- function(calModel,
                         recData,
                         iter = 1000,
                         mixed = FALSE,
                         postcalsamples = NULL,
                         MC = TRUE) {

  if(MC){
    options(mc.cores = parallel::detectCores())
  }else{
    options(mc.cores = 1)
  }

  vects.params <- extract(calModel)

  # Slopes against 10^6/T^2 are positive
  if (mean(vects.params$beta) <= 0) {
    stop("The calibration must use Temperature in 10^6/T^2 (T in K).")
  }

  # Define models
Model <- "
data {
  int<lower=0> n;
  vector[n] y_mes;
  vector[n] y_err;
  int<lower=0> posts;
  vector[posts] alpha;
  vector[posts] beta;
  vector[posts] sigma;
}

parameters {
  matrix<lower=0>[n, posts] x_new;
}

model {
  vector[posts] y_new_hat;
  for(i in 1:n){
    y_new_hat = alpha + beta .* x_new[i,]';
    y_mes[i] ~ normal(y_new_hat, sqrt(sigma^2 + y_err[i]^2));
  }
}

"

totsamp <- if(!mixed){
  length(vects.params[[1]])}else{
    nrow(vects.params[[1]])
}

seqSamples <- if(!is.null(postcalsamples) ){
  sample(1:totsamp, postcalsamples, replace = TRUE)
  }else{
    1:totsamp
  }

# Summarize temperatures (in degrees C) across all posterior draws
sumTemps <- function(fit) {
  Xouts2 <- extract(fit)$x_new
  temps <- sqrt(10 ^ 6 / Xouts2) - 273.15
  t(apply(temps, 2, function(x) c(mean(x), sd(x))))
}

# If mixed perform analysis per material

if (mixed) {
  # Match materials to the labels stored by cal.bayesian()
  mats <- attr(calModel, "Material")
  if (is.null(mats) || anyNA(match(as.character(recData$Material), mats))) {
    stop("Materials in recData must match those in the calibration fitted with cal.bayesian().")
  }
  partMat <- split(recData, recData$Material)

  recs <- lapply(partMat, function(x) {
    j <- match(as.character(x$Material[1]), mats)
    stan_date <- list(
      n = nrow(x),
      y_mes = x$D47,
      y_err = x$D47error,
      posts = length(seqSamples),
      alpha = vects.params$alpha[seqSamples, j],
      beta = vects.params$beta[seqSamples, j],
      sigma = vects.params$sigma[seqSamples]
    )

    data.rstan <- stan(
      data = stan_date,
      model_code = Model,
      chains = 2,
      iter = iter,
      warmup = floor(iter / 2),
      control = list(adapt_delta = 0.90, max_treedepth = 10)
    )

    recs <- sumTemps(data.rstan)

    cbind.data.frame(
      Sample = x$Sample,
      D47 = x$D47,
      D47error = x$D47error,
      meanTemp = recs[, 1],
      error = recs[, 2]
    )
  })
  recs <- do.call(rbind, recs)
  recs <- recs[match(recData$Sample, recs$Sample),]
  row.names(recs) <- NULL
  recs

} else {
  stan_date <- list(
    n = nrow(recData),
    y_mes = recData$D47,
    y_err = recData$D47error,
    posts = length(seqSamples),
    alpha = vects.params$alpha[seqSamples],
    beta = vects.params$beta[seqSamples],
    sigma = vects.params$sigma[seqSamples]
  )

  data.rstan <- stan(
    data = stan_date,
    model_code = Model,
    chains = 2,
    iter = iter,
    warmup = floor(iter / 2)
  )

  recs <- sumTemps(data.rstan)

  reconstructionTable <-
  cbind.data.frame(
    Sample = recData$Sample,
    D47 = recData$D47,
    D47error = recData$D47error,
    meanTemp = recs[, 1],
    error = recs[, 2]
  )
  return(reconstructionTable)
}
}

