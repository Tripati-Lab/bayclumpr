#' Standard error
#' @param x Numeric vector

se <- function(x, iter){
  sqrt(var(x) / iter)
}
