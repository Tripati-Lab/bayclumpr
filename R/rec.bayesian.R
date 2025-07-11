#' This function generate temperature predictions (in 10^6/T2) based on a
#' calibration dataset and target D47. Note that this approach additionally
#' accounts for measured error in the target D47. This approach is congruent
#' with the one used in McClelland et al. (2022).
#'
#' @param calModel The stan model to be analyzed.
#' @param recData The reconstruction dataset.
#' @param iter Number of replicates to retain.
#' @param mixed whether the model \code{calModel} is mixed or not.
#' @param postcalsamples Number of posterior samples to analyze from the calibration step.
#' @param MC Multicore (TRUE/FALSE)
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

  ## Define models
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
  matrix[n, posts] x_new;
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
    1:length(vects.params[[1]])
  }

##If mixed perform analysis per material

if (mixed) {
  partMat <- split(recData, recData$Material)

  recs <- lapply(partMat, function(x) {
    stan_date <- list(
      n = nrow(x),
      y_mes = x$D47,
      y_err = x$D47error,
      posts = length(seqSamples),
      alpha = vects.params$alpha[seqSamples, x$Material[1]],
      beta = vects.params$beta[seqSamples, x$Material[1]],
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

    params2 <- extract(data.rstan)
    Xouts2 <- params2$x_new
    Xdims2 <- dim(Xouts2)
    recs <- sapply(1:iter, function(x){
      rowMeans(Xouts2[x,,])
    })
    recs <- sqrt(10 ^ 6 / recs) - 273.15
    recs <- apply(recs, 1, function(x){
      data.frame(mean(x), sd(x))}
    )

    recs <- do.call(rbind, recs)

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

  params2 <- extract(data.rstan)
  Xouts2 <- params2$x_new
  Xdims2 <- dim(Xouts2)

  recs <- sapply(1:iter, function(x){
    rowMeans(Xouts2[x,,])
  })
  recs <- sqrt(10 ^ 6 / recs) - 273.15
  recs <- apply(recs, 1, function(x){
    data.frame(mean(x), sd(x))}
  )

  recs <- do.call(rbind, recs)

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

