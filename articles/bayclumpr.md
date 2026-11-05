# 1. Getting started

`bayclumpr` fits calibrations between clumped isotopes (D47) and
temperature. It then uses these calibrations to reconstruct temperatures
for new samples. You can fit calibrations under frequentist and Bayesian
frameworks. The same functions run
[BayClump](https://bayclump.tripatilab.epss.ucla.edu/), a Shiny app that
you can use from your browser.

This vignette covers the data you need and a first example. The other
vignettes go into more detail:

1.  Getting started (this vignette).
2.  [`vignette("frequentist", package = "bayclumpr")`](https://tripati-lab.github.io/bayclumpr/articles/frequentist.md):
    OLS, weighted OLS, Deming, and York regressions.
3.  [`vignette("bayesian", package = "bayclumpr")`](https://tripati-lab.github.io/bayclumpr/articles/bayesian.md):
    Bayesian calibrations, priors, model comparison, and
    reconstructions.
4.  [`vignette("covariates", package = "bayclumpr")`](https://tripati-lab.github.io/bayclumpr/articles/covariates.md):
    Bayesian calibrations with a second predictor.

``` r

library(bayclumpr)
```

## The workflow

Every analysis has two steps. First, we fit a calibration. The result is
a table with many plausible values of the intercept (`alpha`) and slope
(`beta`). Second, we use that table to reconstruct temperatures for new
samples. Frequentist functions get these values by bootstrapping.
Bayesian functions get them from the posterior distribution.

| Step | Frequentist | Bayesian |
|:---|:---|:---|
| 1\. Calibrate | [`cal.ols()`](https://tripati-lab.github.io/bayclumpr/reference/cal.ols.md), [`cal.wols()`](https://tripati-lab.github.io/bayclumpr/reference/cal.wols.md), [`cal.deming()`](https://tripati-lab.github.io/bayclumpr/reference/cal.deming.md), [`cal.york()`](https://tripati-lab.github.io/bayclumpr/reference/cal.york.md) | [`cal.bayesian()`](https://tripati-lab.github.io/bayclumpr/reference/cal.bayesian.md), [`cal.ion.bayesian()`](https://tripati-lab.github.io/bayclumpr/reference/cal.ion.bayesian.md) |
| 2\. Get the table of `alpha` and `beta` | Returned by the functions above | [`cal.draws()`](https://tripati-lab.github.io/bayclumpr/reference/cal.draws.md) |
| 3\. Reconstruct | [`rec.clumped()`](https://tripati-lab.github.io/bayclumpr/reference/rec.clumped.md) | [`rec.bayesian()`](https://tripati-lab.github.io/bayclumpr/reference/rec.bayesian.md), [`rec.ion.bayesian()`](https://tripati-lab.github.io/bayclumpr/reference/rec.ion.bayesian.md) |

## Preparing your data

### Calibration data

Calibration datasets need the following columns:

| Column        | Description                                            |
|:--------------|:-------------------------------------------------------|
| `Temperature` | Temperature as 10⁶/T², with T in kelvin                |
| `TempError`   | Uncertainty in `Temperature` (1 SD, in the same units) |
| `D47`         | Measured D47                                           |
| `D47error`    | Uncertainty in `D47` (1 SD)                            |
| `Material`    | Material of each sample. Only needed for mixed models  |

All uncertainties are 1 SD of the value in the table. If `D47` is the
mean of several replicates, `D47error` should be the standard error of
that mean. The same applies to `TempError`.

Temperatures are often reported in °C. To convert them, we use T in
kelvin. For the error, the SD of 10⁶/T² is about 2 × 10⁶ × SD(T) / T³:

``` r

tempC <- c(5, 25, 50)
tempErrorC <- c(0.5, 0.5, 1)
tempK <- tempC + 273.15

data.frame(
  Temperature = 10^6 / tempK^2,
  TempError = 2 * 10^6 * tempErrorC / tempK^3
)
#>   Temperature  TempError
#> 1   12.925337 0.04646894
#> 2   11.249426 0.03773076
#> 3    9.576166 0.05926762
```

### Reconstruction data

Reconstruction datasets need `Sample`, `D47`, and `D47error`. Mixed
models also need `Material`, with the same labels used in the
calibration. Reconstructed temperatures are always returned in °C.

## Simulated data

[`cal.dataset()`](https://tripati-lab.github.io/bayclumpr/reference/cal.dataset.md)
simulates calibration datasets with known parameters (`alpha` = 0.268
and `beta` = 0.0369). There are three error scenarios: low (`"S1"`),
intermediate (`"S2"`), and high (`"S3"`). The columns `x_TRUE` and
`y_TRUE` hold the true values before measurement error.

``` r

set.seed(3)
ds <- cal.dataset(error = "S1", nobs = 50)
head(ds)
#>      x_TRUE Temperature TempError    y_TRUE D47error       D47 Material
#> 1 10.076133   10.089943     0.019 0.6475262   0.0025 0.6442978        1
#> 2 11.414949   11.399569     0.019 0.6841481   0.0025 0.6907357        1
#> 3 12.517576   12.522651     0.019 0.7430624   0.0025 0.7442800        1
#> 4  9.695736    9.662728     0.019 0.6333012   0.0025 0.6354360        1
#> 5 12.391566   12.364749     0.019 0.7379670   0.0025 0.7406881        1
#> 6 12.060248   12.051630     0.019 0.7206252   0.0025 0.7211903        1
```

## A first example

We fit an OLS calibration with 1,000 bootstrap replicates. Each row of
the result is one replicate.

``` r

ols <- cal.ols(data = ds, replicates = 1000)
head(ols)
#>       alpha       beta
#> 1 0.2551112 0.03805364
#> 2 0.2645517 0.03720573
#> 3 0.2739421 0.03638998
#> 4 0.2675038 0.03681503
#> 5 0.2896745 0.03498087
#> 6 0.2682990 0.03664613
colMeans(ols)
#>      alpha       beta 
#> 0.27022607 0.03660694
```

Next, we create a reconstruction dataset for three samples and
reconstruct their temperatures:

``` r

recData <- data.frame(
  Sample = c("Sample 1", "Sample 2", "Sample 3"),
  D47 = c(0.6, 0.7, 0.8),
  D47error = 0.01
)

rec.clumped(recData = recData, obCal = ols)
#>     Sample D47 D47error  meanTemp    error
#> 1 Sample 1 0.6     0.01  60.02587 5.278902
#> 2 Sample 2 0.7     0.01  18.70133 3.477849
#> 3 Sample 3 0.8     0.01 -10.28279 2.663233
```

`meanTemp` is the reconstructed temperature in °C. `error` is its
uncertainty (1 SD). It combines the error in `D47` with the uncertainty
in the calibration.

From here,
[`vignette("frequentist", package = "bayclumpr")`](https://tripati-lab.github.io/bayclumpr/articles/frequentist.md)
compares the four frequentist models.
[`vignette("bayesian", package = "bayclumpr")`](https://tripati-lab.github.io/bayclumpr/articles/bayesian.md)
shows the same workflow under a Bayesian framework.
