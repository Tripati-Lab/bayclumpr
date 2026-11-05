# Temperature reconstructions from bootstrap calibrations

Reconstructs temperatures (in degrees C) from the bootstrap replicates
of a calibration on 10^6/T^2 (T in K) and the measured D47 of each
sample.

## Usage

``` r
rec.clumped(recData, obCal, ndraws = 10000, prediction = FALSE, level = 0.95)
```

## Arguments

- recData:

  Reconstruction dataset. It needs the columns `Sample`, `D47`, and
  `D47error` (1 SD of the D47 value). If D47 is the mean of several
  replicates, `D47error` should be its standard error. Add `N` (number
  of replicates) or `D47sd` (SD of the replicates) to also get `SD`. The
  output of
  [`avg.replicates()`](https://tripati-lab.github.io/bayclumpr/reference/avg.replicates.md)
  has both.

- obCal:

  A `data.frame` with the slopes and intercepts of the calibration, such
  as a result from
  [`cal.ols()`](https://tripati-lab.github.io/bayclumpr/reference/cal.ols.md).
  It also needs `sigma` when `prediction = TRUE`.

- ndraws:

  Number of Monte Carlo draws used for each sample.

- prediction:

  Whether to include the residual scatter of the calibration (`sigma`).
  Use `TRUE` to get the same kind of uncertainty as
  [`rec.bayesian()`](https://tripati-lab.github.io/bayclumpr/reference/rec.bayesian.md).

- level:

  Probability covered by the interval.

## Value

A `data.frame` with one row per sample. `meanTemp` is the temperature
from the mean slope and intercept (in degrees C). `SE` is the standard
error of the temperature, propagated from `D47error`. `SD` is the SD of
the temperature of a single replicate, propagated from the SD of the
replicates (`D47error` x sqrt(`N`), or `D47sd` when `N` is missing). It
is `NA` without replicate information. Both include the uncertainty in
the calibration. `median`, `lower`, and `upper` are the median and the
limits of the interval, based on `SE`.

## Details

The uncertainty always includes `D47error` and the joint spread of
slopes and intercepts in `obCal`. With `prediction = FALSE`, it
describes the temperature predicted by the calibration line. With
`prediction = TRUE`, it also includes the scatter of individual samples
around that line. Draws with 10^6/T^2 at or below zero are not
physically possible. They are removed with a warning.

## Examples

``` r
ds <- cal.dataset(error = "S1", nobs = 50, seed = 1)
ols <- cal.ols(ds, replicates = 100)
recData <- data.frame(Sample = c("A", "B"), D47 = c(0.6, 0.7), D47error = 0.01, N = 3)
rec.clumped(recData, ols)
#>   Sample D47 D47error meanTemp       SD       SE   median    lower    upper
#> 1      A 0.6     0.01 59.87653 9.122010 5.294428 59.84306 50.07710 70.94595
#> 2      B 0.7     0.01 18.61770 5.822095 3.461574 18.67820 12.17715 25.72509
rec.clumped(recData, ols, prediction = TRUE)
#>   Sample D47 D47error meanTemp        SD       SE   median     lower    upper
#> 1      A 0.6     0.01 59.87653 10.882046 8.117084 59.90510 45.500066 77.24463
#> 2      B 0.7     0.01 18.61770  7.255958 5.309079 18.44501  8.764888 29.84227
```
