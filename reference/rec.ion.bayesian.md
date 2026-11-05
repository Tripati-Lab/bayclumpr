# Bayesian temperature reconstructions with a second predictor

Reconstructs temperatures (in degrees C) from the posterior draws of
[`cal.ion.bayesian()`](https://tripati-lab.github.io/bayclumpr/reference/cal.ion.bayesian.md),
the measured D47, and the `Ion` value of each sample.

## Usage

``` r
rec.ion.bayesian(
  calModel,
  recData,
  ndraws = 10000,
  useIonError = FALSE,
  level = 0.95,
  ...
)
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
  `D47error` (1 SD of the D47 value), and `Ion`. Add `N` (number of
  replicates) or `D47sd` (SD of the replicates) to also get `SD`. The
  output of
  [`avg.replicates()`](https://tripati-lab.github.io/bayclumpr/reference/avg.replicates.md)
  has both. It also needs `IonError` (1 SD of `Ion`) when
  `useIonError = TRUE`.

- ndraws:

  Number of draws used for each sample.

- useIonError:

  Whether to account for the error in `Ion`.

- level:

  Probability covered by the interval.

- ...:

  Arguments from earlier versions (`iter`, `postcalsamples`, and `MC`).
  They are ignored.

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

With flat priors on 10^6/T^2 and `Ion`, the posterior of `Ion` follows
its measurement, and the posterior of 10^6/T^2 for each draw is normal.
Both are sampled directly, so no MCMC run is needed. Draws with 10^6/T^2
at or below zero are not physically possible. They are removed with a
warning.

The uncertainty includes the calibration parameters, the residual SD of
the calibration (`sigma`), `D47error`, and `IonError` when used.

## Examples

``` r
# \donttest{
ds <- cal.dataset(error = "S1", nobs = 50, seed = 1)
ds$Ion <- rnorm(nrow(ds))
fit <- cal.ion.bayesian(ds, IonError = 0.1, numSavedSteps = 1000)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'sampling': Boost not found; call install.packages('BH')
recData <- data.frame(Sample = "A", D47 = 0.6, D47error = 0.01, Ion = 0, IonError = 0.1)
rec.ion.bayesian(fit, recData, useIonError = TRUE)
#> Error: object 'fit' not found
# }
```
