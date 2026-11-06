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
cal.draws(fit, model = NULL, newMaterial = FALSE)
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
  `"noerrors"`, `"mixed_errors"`, or `"mixed_noerrors"`, or the name of
  the element (such as `"BLM1_fit"`).

- newMaterial:

  For mixed models, whether to return draws for a material that is not
  in the calibration. These draws combine the uncertainty in the shared
  means with the spread among materials.

## Value

A `data.frame` with one row per draw and the columns `alpha`, `beta`,
and `sigma`. Models with a second predictor also have `gamma`. Mixed
models have one set of draws per material, stacked, with a `Material`
column. With `newMaterial = TRUE`, there is a single set of draws and no
`Material` column.

## Examples

``` r
# \donttest{
ds <- cal.dataset(error = "S1", nobs = 50, seed = 1)
fit <- cal.bayesian(ds, models = "errors", numSavedSteps = 1000)
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
head(cal.draws(fit))
#>       alpha       beta      sigma
#> 1 0.2700057 0.03667362 0.01019158
#> 2 0.2799087 0.03587386 0.01068777
#> 3 0.2660060 0.03657070 0.01263002
#> 4 0.2619720 0.03721020 0.01096346
#> 5 0.2426378 0.03872507 0.01203537
#> 6 0.2849338 0.03548324 0.01145433
# }
```
