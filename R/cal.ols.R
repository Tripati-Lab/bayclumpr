#' Fit OLS regression models on a given calibration dataset
#'
#' @param data The calibration dataset
#' @param replicates Number of bootstrap replicates
#'
#' @importFrom stats lm
#'
#' @export

cal.ols <- function(data, replicates) {
  reps <- lapply(1:replicates, function(x) {
    dataSub <- data[sample(seq_along(data[, 1]), nrow(data), replace = TRUE), ]
    Reg <- summary(lm(D47 ~ Temperature, dataSub))
    res <- cbind.data.frame("alpha" = Reg$coefficients[1, 1], "beta" = Reg$coefficients[2, 1])
    return(res)
  })
  reps <- do.call(rbind, reps)
  return(reps)
}
