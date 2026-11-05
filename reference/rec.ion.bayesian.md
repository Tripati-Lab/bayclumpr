# Bayesian temperature reconstructions with a second predictor

Reconstructs temperatures (in degrees C) from the posterior draws of
[`cal.ion.bayesian()`](https://tripati-lab.github.io/bayclumpr/reference/cal.ion.bayesian.md),
the measured D47, and the `Ion` value of each sample.

## Usage

``` r
rec.ion.bayesian(calModel, recData, ndraws = 10000, useIonError = FALSE, ...)
```

## Arguments

- calModel:

  Posterior draws from
  [`cal.draws()`](https://tripati-lab.github.io/bayclumpr/reference/cal.draws.md),
  a result from
  [`cal.ion.bayesian()`](https://tripati-lab.github.io/bayclumpr/reference/cal.ion.bayesian.md),
  or the `stanfit` in it.

- recData:

  The reconstruction dataset. It needs the columns `Sample`, `D47`,
  `D47error` (1 SD of the D47 value), and `Ion`. It also needs
  `IonError` (1 SD of `Ion`) when `useIonError = TRUE`.

- ndraws:

  Number of draws used for each sample.

- useIonError:

  Whether to account for the error in `Ion`.

- ...:

  Arguments from earlier versions (`iter`, `postcalsamples`, and `MC`).
  They are ignored.

## Value

A `data.frame` with `meanTemp` (posterior mean, in degrees C) and
`error` (1 SD) for each sample. The error includes uncertainty in the
calibration parameters, the residual SD of the calibration (`sigma`),
`D47error`, and `IonError` when used.

## Details

With flat priors on 10^6/T^2 and `Ion`, the posterior of `Ion` follows
its measurement, and the posterior of 10^6/T^2 for each draw is normal.
Both are sampled directly, so no MCMC run is needed. Draws with 10^6/T^2
at or below zero are not physically possible. They are removed with a
warning.
