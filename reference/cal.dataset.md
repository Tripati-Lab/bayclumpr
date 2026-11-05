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
