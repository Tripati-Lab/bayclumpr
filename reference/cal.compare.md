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

In the errors model, the pointwise log-likelihood is conditional on the
latent true temperatures. A warning is given when any Pareto k value is
above 0.7, since the comparison may then be unreliable.
