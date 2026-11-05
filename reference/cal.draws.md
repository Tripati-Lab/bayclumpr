# Posterior draws from a Bayesian calibration

Returns the posterior draws of a Bayesian calibration as a `data.frame`,
in the same shape as the results of the bootstrap calibrations
([`cal.ols()`](https://tripati-lab.github.io/bayclumpr/reference/cal.ols.md)
and others). The result can be passed to
[`rec.bayesian()`](https://tripati-lab.github.io/bayclumpr/reference/rec.bayesian.md),
[`rec.clumped()`](https://tripati-lab.github.io/bayclumpr/reference/rec.clumped.md),
or
[`cal.ci()`](https://tripati-lab.github.io/bayclumpr/reference/cal.ci.md).

## Usage

``` r
cal.draws(fit, model = NULL)
```

## Arguments

- fit:

  A result from
  [`cal.bayesian()`](https://tripati-lab.github.io/bayclumpr/reference/cal.bayesian.md)
  or
  [`cal.ion.bayesian()`](https://tripati-lab.github.io/bayclumpr/reference/cal.ion.bayesian.md),
  or a single `stanfit` from it.

- model:

  The model to use when `fit` holds more than one. Use `"errors"`,
  `"noerrors"`, or `"mixed"`, or the name of the element (such as
  `"BLM1_fit"`).

## Value

A `data.frame` with one row per draw and the columns `alpha`, `beta`,
and `sigma`. Models with a second predictor also have `gamma`. Mixed
models have one set of draws per material, stacked, with a `Material`
column.
