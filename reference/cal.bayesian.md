# Bayesian regressions to calibrate the clumped isotopes paleothermometer using `stan`.

Bayesian regressions to calibrate the clumped isotopes paleothermometer
using `stan`.

## Usage

``` r
cal.bayesian(
  calibrationData,
  numSavedSteps = 3000,
  priors = "Informative",
  MC = TRUE
)
```

## Arguments

- calibrationData:

  The target calibration dataset.

- numSavedSteps:

  Number of MCMC iterations to save.

- priors:

  Either `Informative`, `Weak`, or `Uninformative` on the slope and
  intercept.

- MC:

  Multicore (TRUE/FALSE)
