#' Fit York regression models on a given calibration dataset
#'
#' @param data The calibration dataset
#' @param replicates Number of bootstrap replicates
#' @param samples Number of samples per bootstrap replicate
#'
#' @importFrom IsoplotR york
#'
#' @export

cal.york <- function(data, replicates, samples = NULL) {
  if(is.null(samples)){samples = nrow(data)}
  reps <-
    do.call(rbind, lapply(1:replicates, function(x) {
      dataSub <- data[sample(seq_along(data[, 1]), nrow(data), replace = TRUE), ]
      dataSub$y_SE <- dataSub[, "D47error"]
      dataSub$x_SE <- abs(dataSub$TempError)
      Reg <- york(cbind.data.frame(dataSub$Temperature, dataSub$x_SE, dataSub$D47, dataSub$y_SE))
      cbind.data.frame("alpha" = Reg$a[1], "beta" = Reg$b[1])
    }))
  row.names(reps) <- NULL
  return(reps)
}
