# Summarize a calibration

Summarizes the intercept, slope, and residual SD of a calibration. It
works the same way for bootstrap calibrations
([`cal.ols()`](https://tripati-lab.github.io/bayclumpr/reference/cal.ols.md)
and others) and Bayesian calibrations
([`cal.bayesian()`](https://tripati-lab.github.io/bayclumpr/reference/cal.bayesian.md)
and
[`cal.ion.bayesian()`](https://tripati-lab.github.io/bayclumpr/reference/cal.ion.bayesian.md)).

## Usage

``` r
cal.summary(fit, model = NULL, level = 0.95)
```

## Arguments

- fit:

  A `data.frame` of bootstrap replicates or posterior draws, a result
  from
  [`cal.bayesian()`](https://tripati-lab.github.io/bayclumpr/reference/cal.bayesian.md)
  or
  [`cal.ion.bayesian()`](https://tripati-lab.github.io/bayclumpr/reference/cal.ion.bayesian.md),
  or a single `stanfit` from it.

- model:

  The model to use when `fit` holds more than one, such as `"errors"`.

- level:

  Probability covered by the interval.

## Value

A `data.frame` with one row per parameter and the columns `parameter`,
`estimate` (mean), `SE`, `median`, `lower`, `upper`, `ESS`, and `Rhat`.
Mixed models also have a `Material` column.

## Details

`SE` is the standard error of each estimate. It is the SD of the
bootstrap replicates or of the posterior draws. `sigma` is itself an SD:
the scatter of individual samples around the calibration line. Its `SE`
describes how well that scatter is known.

For Bayesian fits, `ESS` is the effective sample size (`n_eff`) and
`Rhat` the potential scale reduction factor of each parameter, as
reported by `rstan`. Look for `Rhat` at or below 1.01 and at least 100
effective draws per chain. They are `NA` for bootstrap replicates, which
are independent draws, and for a `data.frame` of posterior draws, which
no longer records the chains.

## Examples

``` r
ds <- cal.dataset(error = "S1", nobs = 50, seed = 1)
ols <- cal.ols(ds, replicates = 200)
cal.summary(ols)
#>   parameter   estimate           SE     median       lower      upper ESS Rhat
#> 1     alpha 0.26952294 0.0112427446 0.26848511 0.249811220 0.29615864  NA   NA
#> 2      beta 0.03664000 0.0008750491 0.03672634 0.034499765 0.03827147  NA   NA
#> 3     sigma 0.01202395 0.0013887962 0.01205333 0.009551602 0.01495939  NA   NA
```
