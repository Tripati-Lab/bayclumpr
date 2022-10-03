#' Generate a synthetic dataset for clumped isotopes calibrations
#' @param error Error scenario: low (S1), Intermediate (S2), or High (S3)
#' @param nobs Number of observations in the simulated dataset
#'
#' @importFrom stats rnorm
#'
#' @export

cal.dataset <- function(error = "S1", nobs = 1000){
  set.seed(3)

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

  truex <- rnorm(nobs,12.02585352, 2.5)
  errx <- rnorm(nobs, 0, TempErrC)
  obsx <- truex + errx

  beta1 <- 0.268
  beta2 <- 0.0369
  sdy <- TrueD47
  sdobsy <- AddD

  erry <- rnorm(nobs, 0, sdy)
  truey <- rnorm(nobs,beta1 + beta2*truex,sdobsy)
  obsy <- truey + erry

  ds <- cbind.data.frame(
    x_TRUE = truex,
    Temperature = obsx,
    TempError = errx,
    y_TRUE = truey,
    D47error = erry,
    D47 = obsy,
    Material = 1)

  return(ds)
}


