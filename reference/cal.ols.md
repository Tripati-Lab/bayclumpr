# Fit OLS regression models on a given calibration dataset

Fits an ordinary least squares regression of D47 on temperature to
bootstrap replicates of the calibration dataset. Measurement errors are
not used.

## Usage

``` r
cal.ols(data, replicates = 1000, samples = NULL)
```

## Arguments

- data:

  The calibration dataset. It needs the columns `Temperature` (10^6/T^2,
  with T in K) and `D47`.

- replicates:

  Number of bootstrap replicates

- samples:

  Number of samples per bootstrap replicate. Defaults to the number of
  rows in `data`.

## Value

A `data.frame` with one row per replicate and the columns `alpha`
(intercept), `beta` (slope), and `sigma` (residual SD of D47 around the
fitted line, including measurement error). Failed replicates are removed
with a warning.

## Examples

``` r
ds <- cal.dataset(error = "S1", nobs = 50, seed = 1)
ols <- cal.ols(ds, replicates = 100)
colMeans(ols)
#>      alpha       beta      sigma 
#> 0.26976691 0.03662505 0.01194744 
```
