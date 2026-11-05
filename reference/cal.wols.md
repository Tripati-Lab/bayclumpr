# Fit weighted OLS regression models on a given calibration dataset

Fit weighted OLS regression models on a given calibration dataset

## Usage

``` r
cal.wols(data, replicates, samples = NULL)
```

## Arguments

- data:

  The calibration dataset. It needs the columns `Temperature` (10^6/T^2,
  with T in K) and `D47`.

- replicates:

  Number of bootstrap replicates

- samples:

  Number of samples per bootstrap replicate
