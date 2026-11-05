# Fit weighted OLS regression models on a given calibration dataset

Fits a weighted least squares regression of D47 on temperature to
bootstrap replicates of the calibration dataset. The weights don't come
from the reported errors. They come from a regression of the absolute
OLS residuals on the fitted values, so observations in regions with more
scatter get less weight.

## Usage

``` r
cal.wols(data, replicates = 1000, samples = NULL)
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

## Details

A replicate fails when the regression of the residuals predicts a
non-positive value, since the weights can't be computed.

## Examples

``` r
ds <- cal.dataset(error = "S1", nobs = 50, seed = 1)
wols <- cal.wols(ds, replicates = 100)
colMeans(wols)
#>      alpha       beta      sigma 
#> 0.27136476 0.03649693 0.01195749 
```
