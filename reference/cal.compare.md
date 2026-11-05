# Compare Bayesian calibration models

Compares the models in a
[`cal.bayesian()`](https://tripati-lab.github.io/bayclumpr/reference/cal.bayesian.md)
result using approximate leave-one-out cross-validation (PSIS-LOO).

## Usage

``` r
cal.compare(fit)
```

## Arguments

- fit:

  A result from
  [`cal.bayesian()`](https://tripati-lab.github.io/bayclumpr/reference/cal.bayesian.md)
  with at least two models.

## Value

The result of
[`loo::loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html).
The best model is in the first row.

## Details

In models with measurement errors, the pointwise log-likelihood is
conditional on the latent true temperatures. A warning is given when any
Pareto k value is above 0.7, since the comparison may then be
unreliable.

## Examples

``` r
# \donttest{
ds <- cal.dataset(error = "S1", nobs = 50, seed = 1)
fit <- cal.bayesian(ds, models = c("errors", "noerrors"), numSavedSteps = 1000)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'sampling': Boost not found; call install.packages('BH')
cal.compare(fit)
#> Error: object 'fit' not found
# }
```
