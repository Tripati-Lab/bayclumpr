# Bayesian regressions to calibrate the clumped isotopes paleothermometer using `stan`.

Bayesian regressions to calibrate the clumped isotopes paleothermometer
using `stan`.

## Usage

``` r
cal.ion.bayesian(
  calibrationData,
  numSavedSteps = 3000,
  MC = TRUE,
  IonError = NULL
)
```

## Arguments

- calibrationData:

  The target calibration dataset.

- numSavedSteps:

  Number of MCMC iterations to save.

- MC:

  Multicore (TRUE/FALSE)

- IonError:

  Ion measurement error.
