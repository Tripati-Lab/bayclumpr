#' This function performs temperature reconstructions (in degrees C) from a calibration fit on 10^6/T^2 (T in K).
#'
#' @param recData Reconstruction dataset. It needs the columns \code{Sample}, \code{D47}, and \code{D47error} (1 SD of the D47 value; if D47 is a mean of replicates, use its standard error).
#' @param obCal A \code{data.frame} summarizing the distribution of slopes and intercepts
#' @param ndraws Number of Monte Carlo draws used to estimate the error.
#'
#' @return A \code{data.frame} with \code{meanTemp} (temperature from the mean slope and intercept, in degrees C) and \code{error} (1 SD) for each sample. The error propagates \code{D47error} and the joint spread of slopes and intercepts in \code{obCal}. It does not include the residual scatter of the calibration.
#'
#' @importFrom stats rnorm sd
#'
#' @export

rec.clumped <- function(recData,
                       obCal,
                       ndraws = 10000) {

  # Slopes against 10^6/T^2 are positive
  if (mean(obCal$beta) <= 0) {
    stop("The calibration must use Temperature in 10^6/T^2 (T in K).")
  }

  temp <- sqrt((mean(obCal$beta) * 10^6) /
    (recData$D47 - mean(obCal$alpha))) - 273.15

  # Resample slope and intercept pairs from the calibration
  idx <- sample(nrow(obCal), ndraws, replace = TRUE)
  alpha <- obCal$alpha[idx]
  beta <- obCal$beta[idx]

  # Propagate calibration and D47 uncertainty for each sample
  error <- sapply(seq_len(nrow(recData)), function(i) {
    d47 <- rnorm(ndraws, recData$D47[i], recData$D47error[i])
    sd(sqrt((beta * 10^6) / (d47 - alpha)) - 273.15, na.rm = TRUE)
  })

  recTempS <- cbind.data.frame(
    Sample = recData$Sample,
    D47 = recData$D47,
    D47error = recData$D47error,
    meanTemp = temp,
    error = error
  )
  return(recTempS)
}
