# Bayesian temperature reconstructions

Reconstructs temperatures (in degrees C) from the posterior draws of a
Bayesian calibration and the measured D47 of each sample. This approach
is congruent with the one used in McClelland et al. (2022).

## Usage

``` r
rec.bayesian(calModel, recData, ndraws = 10000, model = NULL, ...)
```

## Arguments

- calModel:

  Posterior draws from
  [`cal.draws()`](https://tripati-lab.github.io/bayclumpr/reference/cal.draws.md),
  a result from
  [`cal.bayesian()`](https://tripati-lab.github.io/bayclumpr/reference/cal.bayesian.md),
  or a single `stanfit` from it.

- recData:

  The reconstruction dataset. It needs the columns `Sample`, `D47`, and
  `D47error` (1 SD of the D47 value). If D47 is the mean of several
  replicates, `D47error` should be its standard error. Mixed models also
  need `Material`, with the same labels used in the calibration.

- ndraws:

  Number of draws used for each sample.

- model:

  The model to use when `calModel` holds more than one, such as
  `"errors"`.

- ...:

  Arguments from earlier versions (`iter`, `mixed`, `postcalsamples`,
  and `MC`). They are ignored.

## Value

A `data.frame` with `meanTemp` (posterior mean, in degrees C) and
`error` (1 SD) for each sample. The error includes uncertainty in the
calibration parameters, the residual SD of the calibration (`sigma`),
and `D47error`. For the errors model, `sigma` excludes the measurement
error of the calibration samples. For the other two models, it includes
it.

## Details

With a flat prior on 10^6/T^2, the posterior of 10^6/T^2 for each
calibration draw is normal. Its mean is (D47 - alpha) / beta and its SD
is sqrt(sigma^2 + D47error^2) / beta. This function samples it directly,
so no MCMC run is needed. Draws with 10^6/T^2 at or below zero are not
physically possible. They are removed with a warning.
