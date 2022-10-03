#' Fit Deming regression models on a given calibration dataset
#'
#' @param data The calibration dataset
#' @param replicates Number of bootstrap replicates
#'
#' @importFrom deming deming
#'
#' @export

cal.deming <- function(data, replicates) {
  reps <-
  do.call(rbind, lapply(1:replicates, function(x) {
    y_SE <- x_SE <- NULL
    dataSub <- data[sample(seq_along(data[, 1]), nrow(data), replace = TRUE), ]
    dataSub$y_SE <- abs(dataSub[, "D47error"])
    dataSub$x_SE <- abs(dataSub$TempError)
    Reg <- deming(D47 ~ Temperature, dataSub, xstd = x_SE, ystd = y_SE)
    cbind.data.frame("alpha" = Reg$coefficients[1], "beta" = Reg$coefficients[2])
  }))
  row.names(reps) <- NULL
  return(reps)
}
