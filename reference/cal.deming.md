# Fit Deming regression models on a given calibration dataset

Fits a Deming regression of D47 on temperature to bootstrap replicates
of the calibration dataset, using
[`deming::deming()`](https://rdrr.io/pkg/deming/man/deming.html). The
errors in both variables are used for each observation.

## Usage

``` r
cal.deming(data, replicates = 1000, samples = NULL)
```

## Arguments

- data:

  The calibration dataset. It needs the columns `Temperature` (10^6/T^2,
  with T in K), `TempError` (1 SD of `Temperature`), `D47`, and
  `D47error` (1 SD of `D47`).

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
dem <- cal.deming(ds, replicates = 100)
colMeans(dem)
#>      alpha       beta      sigma 
#> 0.26851348 0.03672763 0.01194865 
```
