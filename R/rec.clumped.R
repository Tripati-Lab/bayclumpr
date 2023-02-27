#' This function performs temp reconstruction (10^6/T^2 with T in K) for
#' multiple replicates of the same target.
#'
#' @param recData Reconstruction dataset
#' @param obCal A \code{data.frame} summarizing the distribution of slopes and intercepts
#'
#' @export

rec.clumped <- function(recData,
                       obCal) {

  temp <- sqrt((mean(obCal$beta) * 10^6) /
    (recData$D47 - mean(obCal$alpha))) - 273.15

  temp_E <- sqrt((mean(obCal$beta) * 10^6) /
    (recData$D47 + recData$D47error - mean(obCal$alpha))) - 273.15
  error <- (temp - temp_E)

  recTempS <- cbind.data.frame(
    Sample = recData$Sample,
    D47 = recData$D47,
    D47error = recData$D47error,
    meanTemp = temp,
    error = error
  )
  return(recTempS)
}
