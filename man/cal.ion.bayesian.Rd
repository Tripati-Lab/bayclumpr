% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cal.ion.bayesian.R
\name{cal.ion.bayesian}
\alias{cal.ion.bayesian}
\title{Bayesian regressions to calibrate the clumped isotopes paleothermometer using
\code{stan}.}
\usage{
cal.ion.bayesian(
  calibrationData,
  numSavedSteps = 3000,
  MC = TRUE,
  useIonError = FALSE
)
}
\arguments{
\item{calibrationData}{The target calibration dataset.}

\item{numSavedSteps}{Number of MCMC iterations to save.}

\item{MC}{Multicore (TRUE/FALSE)}

\item{useIonError}{Whether or not to account for ion measurement error.}
}
\description{
Bayesian regressions to calibrate the clumped isotopes paleothermometer using
\code{stan}.
}
