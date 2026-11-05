#' Generate a synthetic dataset for clumped isotopes calibrations
#'
#' Simulates a calibration dataset with alpha = 0.268 and beta = 0.0369. True temperatures (10^6/T^2) are drawn from a normal distribution with mean 12 and SD 2.
#'
#' @param error Error scenario: low (S1), Intermediate (S2), or High (S3)
#' @param nobs Number of observations in the simulated dataset
#' @param seed Optional seed. When given, it is passed to \code{set.seed()}, which also sets the seed for the rest of the session.
#'
#' @return A \code{data.frame} with the true and measured values of \code{Temperature} and \code{D47}, their errors (1 SD), and \code{Material} (always 1).
#'
#' @importFrom stats rnorm
#'
#' @export

cal.dataset <- function(error = "S1", nobs = 1000, seed = NULL){

  error <- match.arg(error, c("S1", "S2", "S3"))
  if (!is.null(seed)) set.seed(seed)

  if(error == 'S1'){
    #Low
    AddD = 0.0025 #Bernasconi
    TrueD47 = 0.0125 #Petersen
    TempErrC = 0.019
  }

  if(error == 'S2'){
    #Intermediate
    AddD = 0.0075
    TrueD47 = 0.0225
    TempErrC = 0.077
  }

  if(error == 'S3'){
    #High
    AddD = 0.0125
    TrueD47 = 0.0275
    TempErrC = 0.155
  }


  # simulate covariate data
  n <- nobs
  sdx <- 2
  pop_taux <- 1 / (sdx * sdx)
  sdobs <- TempErrC
  taux <- 1 / (sdobs * sdobs)
  truex <- rnorm(n, 12, sdx)
  errorx <- rnorm(n, 0, sdobs)
  obsx <- truex + errorx

  # simulate response data
  alpha <- 0.268
  beta <- 0.0369

  # process error
  p_sdy <-  TrueD47
  p_tauy <- 1 / (p_sdy * p_sdy)
  p_errory <- rnorm(n, 0, p_sdy)

  # measurement error
  obs_sdy <- AddD
  obs_tauy <- 1 / (obs_sdy * obs_sdy)
  obs_errory <- rnorm(n,0,obs_sdy)

  truey <- alpha + beta*truex + p_errory
  obsy <- truey + obs_errory

  ds <- cbind.data.frame(
    x_TRUE = truex,
    Temperature = obsx,
    TempError = TempErrC,
    y_TRUE = truey,
    D47error = AddD,
    D47 = obsy,
    Material = 1)

  return(ds)
}


