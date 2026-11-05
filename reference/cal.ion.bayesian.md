# Bayesian calibration with a second predictor

Fits a Bayesian regression of D47 on temperature and a second predictor
(`Ion`) using Stan.

## Usage

``` r
cal.ion.bayesian(
  calibrationData,
  IonError = NULL,
  numSavedSteps = 3000,
  chains = 4,
  cores = getOption("mc.cores", 1L),
  MC = NULL
)
```

## Arguments

- calibrationData:

  The calibration dataset. It needs the columns `Temperature` (10^6/T^2,
  with T in K), `TempError` (1 SD of `Temperature`), `D47`, `D47error`
  (1 SD of `D47`), and `Ion`.

- IonError:

  Ion measurement error (1 SD). Use a single value for all samples or
  one value per sample. Use `NULL` to treat `Ion` as known without
  error.

- numSavedSteps:

  Total number of posterior draws to keep across all chains.

- chains:

  Number of MCMC chains.

- cores:

  Number of cores used to run the chains in parallel.

- MC:

  Deprecated. Use `cores` instead.

## Value

An object of class `bayclump_cal` with a single `stanfit` named
`BLM1_fit`. Use
[`cal.draws()`](https://tripati-lab.github.io/bayclumpr/reference/cal.draws.md)
to get the posterior draws.

## Details

The model accounts for errors in `Temperature` and `D47`, and in `Ion`
when `IonError` is given. Its `sigma` is the scatter left after removing
the measurement error in `D47`. The true values of each predictor have a
normal prior with a mean and SD estimated from the data. The intercept
(alpha) and both slopes (beta and gamma) have flat priors.
