# Bayesian temperature reconstructions

Reconstructs temperatures (in degrees C) from the posterior draws of a
Bayesian calibration and the measured D47 of each sample. This approach
is congruent with the one used in McClelland et al. (2022).

## Usage

``` r
rec.bayesian(
  calModel,
  recData,
  ndraws = 10000,
  model = NULL,
  newMaterial = FALSE,
  level = 0.95,
  ...
)
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
  replicates, `D47error` should be its standard error. Add `N` (number
  of replicates) or `D47sd` (SD of the replicates) to also get `SD`. The
  output of
  [`avg.replicates()`](https://tripati-lab.github.io/bayclumpr/reference/avg.replicates.md)
  has both. Mixed models also need `Material`, with the same labels used
  in the calibration.

- ndraws:

  Number of draws used for each sample.

- model:

  The model to use when `calModel` holds more than one, such as
  `"errors"`.

- newMaterial:

  For mixed models, whether the samples belong to a material that is not
  in the calibration. When `TRUE`, `Material` in `recData` is ignored.
  When `calModel` is a table of draws, get them with
  `cal.draws(fit, newMaterial = TRUE)` instead.

- level:

  Probability covered by the interval.

- ...:

  Arguments from earlier versions (`iter`, `mixed`, `postcalsamples`,
  and `MC`). They are ignored.

## Value

A `data.frame` with one row per sample. `meanTemp` is the posterior mean
temperature (in degrees C). `SE` is the standard error of the
temperature, propagated from `D47error`. `SD` is the SD of the
temperature of a single replicate, propagated from the SD of the
replicates (`D47error` x sqrt(`N`), or `D47sd` when `N` is missing). It
is `NA` without replicate information. Both include the uncertainty in
the calibration. `median`, `lower`, and `upper` are the posterior median
and the limits of the interval, based on `SE`.

## Details

With a flat prior on 10^6/T^2, the posterior of 10^6/T^2 for each
calibration draw is normal. Its mean is (D47 - alpha) / beta and its SD
is sqrt(sigma^2 + D47error^2) / beta. This function samples it directly,
so no MCMC run is needed. Draws with 10^6/T^2 at or below zero are not
physically possible. They are removed with a warning.

The uncertainty includes the calibration parameters, the residual SD of
the calibration (`sigma`), and `D47error`. For models with measurement
errors, `sigma` excludes the measurement error of the calibration
samples. For models without them, it includes it.

## Examples

``` r
# \donttest{
ds <- cal.dataset(error = "S1", nobs = 50, seed = 1)
fit <- cal.bayesian(ds, models = "errors", numSavedSteps = 1000)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'sampling': Boost not found; call install.packages('BH')
recData <- data.frame(Sample = c("A", "B"), D47 = c(0.6, 0.7), D47error = 0.01, N = 3)
rec.bayesian(fit, recData)
#> Error: object 'fit' not found
# }
```
