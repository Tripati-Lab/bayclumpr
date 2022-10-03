#' Standard error
#' @param x Numeric vector
#' @param iter Number of iterations
#'
#' @importFrom stats var

se <- function(x, iter){
  sqrt(var(x, na.rm = TRUE) / iter)
}
