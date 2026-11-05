# Fit York regression models on a given calibration dataset

Fit York regression models on a given calibration dataset

## Usage

``` r
cal.york(data, replicates, samples = NULL)
```

## Arguments

- data:

  The calibration dataset. It needs the columns `Temperature` (10^6/T^2,
  with T in K), `TempError` (1 SD of `Temperature`), `D47`, and
  `D47error` (1 SD of `D47`).

- replicates:

  Number of bootstrap replicates

- samples:

  Number of samples per bootstrap replicate
