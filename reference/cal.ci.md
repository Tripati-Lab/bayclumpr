# Calibration line with its interval

Computes the median calibration line and its interval over a range of
temperatures. This is mostly used for plotting.

## Usage

``` r
cal.ci(data, from, to, length.out = 100, level = 0.95)
```

## Arguments

- data:

  A `data.frame` with the columns `alpha` and `beta`, such as a result
  from
  [`cal.ols()`](https://tripati-lab.github.io/bayclumpr/reference/cal.ols.md)
  or
  [`cal.draws()`](https://tripati-lab.github.io/bayclumpr/reference/cal.draws.md).
  For mixed models, use the rows of a single material.

- from:

  The lowest temperature (10^6/T^2).

- to:

  The highest temperature (10^6/T^2).

- length.out:

  Number of temperatures between `from` and `to`.

- level:

  Probability covered by the interval.

## Value

A `data.frame` with the columns `x` (temperature), `median_est` (median
D47), and `ci_lower_est` and `ci_upper_est` (limits of the interval).

## Examples

``` r
ds <- cal.dataset(error = "S1", nobs = 50, seed = 1)
ols <- cal.ols(ds, replicates = 100)
ci <- cal.ci(ols, from = min(ds$Temperature), to = max(ds$Temperature))
head(ci)
#>          x median_est ci_lower_est ci_upper_est
#> 1 7.571132  0.5471813    0.5379169    0.5556029
#> 2 7.647879  0.5499803    0.5408654    0.5582916
#> 3 7.724627  0.5527793    0.5438140    0.5609802
#> 4 7.801374  0.5555805    0.5467626    0.5636689
#> 5 7.878121  0.5583910    0.5497112    0.5663575
#> 6 7.954868  0.5612052    0.5526597    0.5690462
```
