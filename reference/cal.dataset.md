# Generate a synthetic dataset for clumped isotopes calibrations

Simulates a calibration dataset with alpha = 0.268 and beta = 0.0369.
True temperatures (10^6/T^2) are drawn from a normal distribution with
mean 12 and SD 2.

## Usage

``` r
cal.dataset(error = "S1", nobs = 1000, seed = NULL)
```

## Arguments

- error:

  Error scenario: low (S1), Intermediate (S2), or High (S3)

- nobs:

  Number of observations in the simulated dataset

- seed:

  Optional seed. When given, it is passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html), which also sets
  the seed for the rest of the session.

## Value

A `data.frame` with the true and measured values of `Temperature` and
`D47`, their errors (1 SD), and `Material` (always 1).

## Examples

``` r
ds <- cal.dataset(error = "S2", nobs = 50, seed = 1)
head(ds)
#>     x_TRUE Temperature TempError    y_TRUE D47error       D47 Material
#> 1 10.74709    10.77775     0.077 0.6506095   0.0075 0.6539859        1
#> 2 12.36729    12.32016     0.077 0.7253005   0.0075 0.7251613        1
#> 3 10.32874    10.35501     0.077 0.6286349   0.0075 0.6262494        1
#> 4 15.19056    15.10360     0.077 0.8320874   0.0075 0.8251172        1
#> 5 12.65902    12.76936     0.077 0.7203895   0.0075 0.7092336        1
#> 6 10.35906    10.51155     0.077 0.6900134   0.0075 0.6819495        1
```
