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
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Compiling the Stan model noerrors. This happens once per session and can take a minute or two.
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> Warning: Check convergence for: BLM1_fit_NoErrors. Print the result for details.
cal.compare(fit)
#>              model elpd_diff se_diff p_worse diag_diff diag_elpd
#>  BLM1_fit_NoErrors       0.0     0.0      NA                    
#>           BLM1_fit       0.0     0.1    0.61   N < 100          
#> 
#> Diagnostic flags present.
#> See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
#> or https://mc-stan.org/loo/reference/loo-glossary.html.
# }
```
