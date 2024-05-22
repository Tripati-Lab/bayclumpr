#' Fit OLS regression models on a given calibration dataset
#'
#' @param data The calibration dataset
#' @param replicates Number of bootstrap replicates
#' @param samples Number of samples per bootstrap replicate
#'
#' @importFrom stats lm
#'
#' @export

cal.ols <- function(data, replicates, samples = NULL) {
  if(is.null(samples)){samples = nrow(data)}
  reps <- lapply(1:replicates, function(x) {
    tryCatch({
    dataSub <- data[sample(seq_along(data[, 1]), nrow(data), replace = TRUE), ]
    Reg <- summary(lm(D47 ~ Temperature, dataSub))
    res <- cbind.data.frame("alpha" = Reg$coefficients[1, 1], "beta" = Reg$coefficients[2, 1])
    return(res)
    }, error=function(e){})
  })
  reps <- do.call(rbind, reps)
  return(reps)
}
