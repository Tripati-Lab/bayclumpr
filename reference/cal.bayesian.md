# Bayesian calibration of the clumped isotope paleothermometer

Fits up to three Bayesian linear regressions of D47 on temperature using
Stan.

## Usage

``` r
cal.bayesian(
  calibrationData,
  models = c("errors", "noerrors", "mixed"),
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
  `D47error` (1 SD of `D47`). The mixed model also needs `Material`.

- models:

  The models to fit. Use any of `"errors"`, `"noerrors"`, and `"mixed"`.

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

  Number of decimal places to print.

- ...:

  Not used.

## Value

An object of class `bayclump_cal`. It is a list of `stanfit` objects
named `BLM1_fit`, `BLM1_fit_NoErrors`, and `BLM3_fit`, with the priors
and a convergence summary stored as attributes. Use
[`cal.draws()`](https://tripati-lab.github.io/bayclumpr/reference/cal.draws.md)
to get the posterior draws and
[`cal.compare()`](https://tripati-lab.github.io/bayclumpr/reference/cal.compare.md)
to compare models.

## Details

The three models are:

- `"errors"` (`BLM1_fit`) accounts for errors in both `Temperature` and
  `D47`. Its `sigma` is the scatter left after removing the measurement
  error in `D47`. The true temperatures have a normal prior with a mean
  and SD estimated from the data.

- `"noerrors"` (`BLM1_fit_NoErrors`) ignores both errors, so its `sigma`
  includes the measurement error.

- `"mixed"` (`BLM3_fit`) lets intercepts and slopes vary by `Material`
  around shared means (partial pooling). It ignores measurement errors.
  The priors apply to the shared means. The SDs among materials have
  half-normal priors with scales of 0.1 (intercepts) and 0.01 (slopes).
  With fewer than five materials, these SDs are mostly set by their
  priors. The model is skipped when there are fewer than two materials.

The `"Informative"` values are the ones used in earlier versions of this
package. Check that they suit the reference frame of your data, or set
your own.

A warning is given when a model shows signs of poor convergence: Rhat
above 1.01, fewer than 100 effective draws per chain, or divergent
transitions.
