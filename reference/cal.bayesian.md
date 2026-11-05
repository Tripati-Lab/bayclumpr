# Bayesian calibration of the clumped isotope paleothermometer

Fits up to four Bayesian linear regressions of D47 on temperature using
Stan. The models cover every combination of two choices: whether
measurement errors are used, and whether intercepts and slopes vary by
material.

## Usage

``` r
cal.bayesian(
  calibrationData,
  models = c("errors", "noerrors", "mixed_errors", "mixed_noerrors"),
  priors = "Informative",
  numSavedSteps = 3000,
  chains = 4,
  cores = getOption("mc.cores", 1L),
  MC = NULL
)

# S3 method for class 'bayclump_cal'
print(x, digits = 4, ...)
```

## Arguments

- calibrationData:

  The calibration dataset. It needs the columns `Temperature` (10^6/T^2,
  with T in K), `TempError` (1 SD of `Temperature`), `D47`, and
  `D47error` (1 SD of `D47`). Mixed models also need `Material`.

- models:

  The models to fit. Use any of `"errors"`, `"noerrors"`,
  `"mixed_errors"`, and `"mixed_noerrors"`. `"mixed"` is kept as another
  name for `"mixed_noerrors"`.

- priors:

  Priors on the intercept (alpha) and slope (beta). Use `"Informative"`,
  `"Weak"`, or `"Flat"`, or a list such as
  `list(alpha = c(0.231, 0.065), beta = c(0.039, 0.004))` with the mean
  and SD of each normal prior. `"Uninformative"` is kept as another name
  for `"Flat"`.

- numSavedSteps:

  Total number of posterior draws to keep across all chains.

- chains:

  Number of MCMC chains.

- cores:

  Number of cores used to run the chains in parallel.

- MC:

  Deprecated. Use `cores` instead.

- x:

  A `bayclump_cal` object.

- digits:

  Number of significant digits to print.

- ...:

  Not used.

## Value

An object of class `bayclump_cal`. It is a list of `stanfit` objects
named as in the table above, with the priors and a convergence summary
stored as attributes. Printing it shows
[`cal.summary()`](https://tripati-lab.github.io/bayclumpr/reference/cal.summary.md)
for each model and the convergence summary. Use
[`cal.summary()`](https://tripati-lab.github.io/bayclumpr/reference/cal.summary.md)
to get these tables,
[`cal.draws()`](https://tripati-lab.github.io/bayclumpr/reference/cal.draws.md)
to get the posterior draws, and
[`cal.compare()`](https://tripati-lab.github.io/bayclumpr/reference/cal.compare.md)
to compare models.

## Details

The four models are:

|  |  |  |  |
|----|----|----|----|
| `models` | Name in the result | Measurement errors | Varies by `Material` |
| `"errors"` | `BLM1_fit` | Yes | No |
| `"noerrors"` | `BLM1_fit_NoErrors` | No | No |
| `"mixed_errors"` | `BLM3_fit_Errors` | Yes | Yes |
| `"mixed_noerrors"` | `BLM3_fit` | No | Yes |

Models with measurement errors use `TempError` and `D47error` for each
observation. Their `sigma` is the scatter left after removing the
measurement error in `D47`. The true temperatures have a normal prior
with a mean and SD estimated from the data. Models without measurement
errors ignore both errors, so their `sigma` includes the measurement
error.

Mixed models let intercepts and slopes vary by `Material` around shared
means (partial pooling). The priors set by `priors` apply to the shared
means. The SDs among materials have half-normal priors with scales of
0.1 for the intercept at the mean calibration temperature and 0.01 for
the slope. With fewer than five materials, these SDs are mostly set by
their priors. Mixed models are skipped when there are fewer than two
materials.

Temperature is centered inside every model to improve sampling. Results
are reported for the usual intercept at 10^6/T^2 = 0, where the priors
are also set.

The `"Informative"` values are the ones used in earlier versions of this
package. Check that they suit the reference frame of your data, or set
your own.

A warning is given when a model shows signs of poor convergence: Rhat
above 1.01, fewer than 100 effective draws per chain, or divergent
transitions.

## Examples

``` r
# \donttest{
ds <- cal.dataset(error = "S1", nobs = 50, seed = 1)
fit <- cal.bayesian(ds, models = "errors", numSavedSteps = 1000)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'sampling': Boost not found; call install.packages('BH')
fit
#> Error: object 'fit' not found
# }
```
