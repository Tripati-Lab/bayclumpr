# This function performs temperature reconstructions (in degrees C) from a calibration fit on 10^6/T^2 (T in K).

This function performs temperature reconstructions (in degrees C) from a
calibration fit on 10^6/T^2 (T in K).

## Usage

``` r
rec.clumped(recData, obCal, ndraws = 10000)
```

## Arguments

- recData:

  Reconstruction dataset. It needs the columns `Sample`, `D47`, and
  `D47error` (1 SD of the D47 value). If D47 is the mean of several
  replicates, `D47error` should be its standard error.

- obCal:

  A `data.frame` summarizing the distribution of slopes and intercepts

- ndraws:

  Number of Monte Carlo draws used to estimate the error.

## Value

A `data.frame` with `meanTemp` (temperature from the mean slope and
intercept, in degrees C) and `error` (1 SD) for each sample. The error
propagates `D47error` and the joint spread of slopes and intercepts in
`obCal`. It does not include the residual scatter of the calibration.
