# 4. Calibrations with a second predictor

Some calibrations depend on a second variable, such as the concentration
of an ion. This vignette shows how to include it as a second predictor:

D47 = `alpha` + `beta` × `Temperature` + `gamma` × `Ion`

It builds on
[`vignette("bayesian", package = "bayclumpr")`](https://tripati-lab.github.io/bayclumpr/articles/bayesian.md).

``` r

library(bayclumpr)
```

## Preparing the data

The calibration dataset needs the usual columns plus `Ion`. If `Ion` is
measured with error, we also need its uncertainty (1 SD). We simulate
100 observations where `Ion` raises D47 by 0.005 per unit. We then add
measurement error to `Ion` with an SD of 0.1.

``` r

set.seed(3)
ds <- cal.dataset(error = "S1", nobs = 100)

trueIon <- rnorm(nrow(ds))
ds$D47 <- ds$D47 + 0.005 * trueIon
ds$IonError <- 0.1
ds$Ion <- trueIon + rnorm(nrow(ds), 0, ds$IonError)
```

## Fitting the calibration

[`cal.ion.bayesian()`](https://tripati-lab.github.io/bayclumpr/reference/cal.ion.bayesian.md)
accounts for errors in `Temperature`, `D47`, and `Ion`. `IonError` can
be a single value for all samples or one value per sample.

``` r

fit <- cal.ion.bayesian(calibrationData = ds, IonError = ds$IonError)
fit
#> Bayesian clumped isotope calibration
#> Priors: Flat 
#> 
#> BLM1_fit
#>  parameter estimate        SE   median    lower    upper
#>      alpha 0.272221 0.0101978 0.272289 0.252790 0.292457
#>       beta 0.036662 0.0008437 0.036654 0.035025 0.038262
#>      gamma 0.005016 0.0012738 0.005045 0.002435 0.007419
#>      sigma 0.013616 0.0010132 0.013565 0.011776 0.015795
#> 
#> Convergence
#>     model max_Rhat min_ESS divergences
#>  BLM1_fit   0.9997    2474           0
```

`gamma` is the slope for `Ion`. Its posterior should be close to the
simulated value of 0.005.

When `Ion` is known without error, we leave out `IonError`:

``` r

fitKnown <- cal.ion.bayesian(calibrationData = ds)
fitKnown
#> Bayesian clumped isotope calibration
#> Priors: Flat 
#> 
#> BLM1_fit
#>  parameter estimate        SE   median    lower    upper
#>      alpha 0.272166 0.0097646 0.272049 0.253372 0.291367
#>       beta 0.036670 0.0008032 0.036674 0.035068 0.038200
#>      gamma 0.005044 0.0012516 0.005059 0.002507 0.007482
#>      sigma 0.013589 0.0010361 0.013542 0.011749 0.015711
#> 
#> Convergence
#>     model max_Rhat min_ESS divergences
#>  BLM1_fit   0.9997    2634           0
```

## Reconstructing temperatures

The reconstruction dataset needs `Ion` for each sample, and `IonError`
when we account for its uncertainty:

``` r

recData <- data.frame(
  Sample = paste("Sample", 1:3),
  D47 = c(0.6, 0.7, 0.8),
  D47error = 0.01,
  Ion = c(-1, 0, 1),
  IonError = 0.1
)
```

With `useIonError = TRUE`, the error in `Ion` is included in the
reconstructed temperatures:

``` r

rec.ion.bayesian(calModel = fit, recData = recData, useIonError = TRUE)
#>     Sample D47 D47error  meanTemp SD       SE   median      lower      upper
#> 1 Sample 1 0.6     0.01 59.156851 NA 8.567143 58.98611  43.048420 76.4799251
#> 2 Sample 2 0.7     0.01 19.920086 NA 5.793169 19.65157   9.238422 32.0191140
#> 3 Sample 3 0.8     0.01 -8.176717 NA 4.342052 -8.29918 -16.420039  0.6563953
```

With `useIonError = FALSE`, `Ion` is treated as known:

``` r

rec.ion.bayesian(calModel = fit, recData = recData, useIonError = FALSE)
#>     Sample D47 D47error  meanTemp SD       SE    median     lower     upper
#> 1 Sample 1 0.6     0.01 59.118286 NA 8.613112 58.910189  43.18238 77.029470
#> 2 Sample 2 0.7     0.01 19.722917 NA 5.792888 19.650943   8.92384 31.400810
#> 3 Sample 3 0.8     0.01 -8.261464 NA 4.338728 -8.299877 -16.58234  0.606063
```

The columns are the same as in
[`rec.bayesian()`](https://tripati-lab.github.io/bayclumpr/reference/rec.bayesian.md).
`meanTemp` is the posterior mean temperature in °C, and `SE` is its
standard error. `SD` is `NA` here because the dataset has neither `N`
nor `D47sd`. Add one of them to get the uncertainty of a single
replicate. `median`, `lower`, and `upper` describe the 95% posterior
interval, based on `SE`. The uncertainty includes the calibration, the
scatter around the calibration line, and the error in `D47`.
