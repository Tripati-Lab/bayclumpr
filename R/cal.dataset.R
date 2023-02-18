#' Generate a synthetic dataset for clumped isotopes calibrations
#' @param error Error scenario: low (S1), Intermediate (S2), or High (S3)
#' @param nobs Number of observations in the simulated dataset
#'
#' @importFrom stats rnorm
#' @importFrom EIVmodels sim_slr
#'
#' @export

cal.dataset <- function(error = "S1", nobs = 1000){
  set.seed(3)

  if(error == 'S1'){
    #Low
    AddD= 0.0025 #Bernasconi
    TrueD47= 0.0125 #Petersen
    TempErrC= 0.019
  }

  if(error == 'S2'){
    #Intermediate
    AddD= 0.0075
    TrueD47= 0.0225
    TempErrC= 0.077
  }

  if(error == 'S3'){
    #High
    AddD= 0.0125
    TrueD47= 0.0275
    TempErrC= 0.155
  }


  data <- sim_slr(n_sim = nobs,
                  alpha = 0.268,
                  beta = 0.0369,
                  y_err = TrueD47,
                  x_err = TempErrC,
                  min_x = 5,
                  max_x = 19)

  ds <- cbind.data.frame(
    x_TRUE = data$true_x,
    Temperature = data$x,
    TempError = data$x_err,
    y_TRUE = data$true_y,
    D47error = data$y_err,
    D47 = data$y,
    Material = 1)

  return(ds)
}


