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
(alpha) and both slopes (beta and gamma) have flat priors. Both
predictors are centered inside the model to improve sampling, and alpha
is reported for `Temperature` = 0 and `Ion` = 0.

## Examples

``` r
# \donttest{
ds <- cal.dataset(error = "S1", nobs = 50, seed = 1)
ds$Ion <- rnorm(nrow(ds))
fit <- cal.ion.bayesian(ds, IonError = 0.1, numSavedSteps = 1000)
#> Compiling the Stan model ion_errors. This happens once per session and can take a minute or two.
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
fit
#> Bayesian clumped isotope calibration
#> Priors: Flat 
#> 
#> BLM1_fit
#>  parameter estimate       SE   median    lower    upper    ESS   Rhat
#>      alpha 0.269171 0.012994 0.269376  0.24379 0.293605  982.9 0.9982
#>       beta 0.036663 0.001057 0.036662  0.03467 0.038730  969.8 0.9987
#>      gamma 0.001084 0.001651 0.001072 -0.00210 0.004249 1147.6 1.0007
#>      sigma 0.012306 0.001364 0.012178  0.01001 0.015404 2290.9 0.9981
#> 
#> Convergence
#>     model max_Rhat min_ESS divergences
#>  BLM1_fit    1.001   969.8           0
# }
```
