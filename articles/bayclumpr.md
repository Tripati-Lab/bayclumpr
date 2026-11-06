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
5.  [`vignette("comparison", package = "bayclumpr")`](https://tripati-lab.github.io/bayclumpr/articles/comparison.md):
    frequentist and Bayesian approaches on the same dataset, with plots
    of parameters and reconstructions.

``` r

library(bayclumpr)
```

## The workflow

Every analysis has two steps. First, we fit a calibration. The result is
a table with many plausible values of the intercept (`alpha`), the slope
(`beta`), and the scatter around the calibration line (`sigma`). Second,
we use that table to reconstruct temperatures for new samples.
Frequentist functions get these values by bootstrapping. Bayesian
functions get them from the posterior distribution.

| Step | Frequentist | Bayesian |
|:---|:---|:---|
| 1\. Calibrate | [`cal.ols()`](https://tripati-lab.github.io/bayclumpr/reference/cal.ols.md), [`cal.wols()`](https://tripati-lab.github.io/bayclumpr/reference/cal.wols.md), [`cal.deming()`](https://tripati-lab.github.io/bayclumpr/reference/cal.deming.md), [`cal.york()`](https://tripati-lab.github.io/bayclumpr/reference/cal.york.md) | [`cal.bayesian()`](https://tripati-lab.github.io/bayclumpr/reference/cal.bayesian.md), [`cal.ion.bayesian()`](https://tripati-lab.github.io/bayclumpr/reference/cal.ion.bayesian.md) |
| 2\. Get the table of `alpha`, `beta`, and `sigma` | Returned by the functions above | [`cal.draws()`](https://tripati-lab.github.io/bayclumpr/reference/cal.draws.md) |
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

### Converting temperatures

Temperatures are often reported in °C.
[`temp2x()`](https://tripati-lab.github.io/bayclumpr/reference/temp2x.md)
converts them to 10⁶/T² and propagates their errors. The result has the
column names that the calibration functions expect:

``` r

temp2x(c(5, 25, 50), tempError = c(0.5, 0.5, 1))
#>   Temperature  TempError
#> 1   12.925337 0.04646894
#> 2   11.249426 0.03773076
#> 3    9.576166 0.05926762
```

[`x2temp()`](https://tripati-lab.github.io/bayclumpr/reference/temp2x.md)
converts them back to °C.

### Averaging replicates

When each sample has several replicate measurements,
[`avg.replicates()`](https://tripati-lab.github.io/bayclumpr/reference/avg.replicates.md)
averages them and returns the standard error of each mean as `D47error`:

``` r

reps <- data.frame(
  Sample = rep(c("A", "B"), each = 3),
  D47 = c(0.601, 0.598, 0.605, 0.702, 0.695, 0.699),
  Material = "calcite"
)
avg.replicates(reps)
#>   Sample       D47       D47sd N    D47error Material
#> 1      A 0.6013333 0.003511885 3 0.002027588  calcite
#> 2      B 0.6986667 0.003511885 3 0.002027588  calcite
```

By default, the SE comes from the SD of each sample’s replicates. With
few replicates, this SD is poorly constrained. The `pooledSD` argument
uses a common SD instead, such as the long-term SD of standards measured
in the same lab.

### Reconstruction data

Reconstruction datasets need `Sample`, `D47`, and `D47error`. Mixed
models also need `Material`, with the same labels used in the
calibration.

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
#>       alpha       beta      sigma
#> 1 0.2551112 0.03805364 0.01270688
#> 2 0.2645517 0.03720573 0.01338969
#> 3 0.2739421 0.03638998 0.01147219
#> 4 0.2675038 0.03681503 0.01205161
#> 5 0.2896745 0.03498087 0.01064595
#> 6 0.2682990 0.03664613 0.01425720
```

[`cal.summary()`](https://tripati-lab.github.io/bayclumpr/reference/cal.summary.md)
summarizes any calibration in the same way:

``` r

cal.summary(ols)
#>   parameter   estimate           SE     median      lower      upper ESS Rhat
#> 1     alpha 0.27022607 0.0106035367 0.27029091 0.24894691 0.29123768  NA   NA
#> 2      beta 0.03660694 0.0009469391 0.03662508 0.03469022 0.03850465  NA   NA
#> 3     sigma 0.01298611 0.0011871156 0.01299158 0.01067024 0.01526967  NA   NA
```

`estimate` is the mean of each parameter and `SE` is its standard error,
which is the SD of the bootstrap replicates. `sigma` is itself an SD:
the scatter of individual samples around the calibration line.

Next, we create a reconstruction dataset for three samples, each
measured three times, and reconstruct their temperatures:

``` r

recData <- data.frame(
  Sample = c("Sample 1", "Sample 2", "Sample 3"),
  D47 = c(0.6, 0.7, 0.8),
  D47error = 0.006,
  N = 3
)

rec.clumped(recData = recData, obCal = ols)
#>     Sample D47 D47error  meanTemp       SD       SE    median     lower
#> 1 Sample 1 0.6    0.006  60.02587 5.451499 3.329676  59.99816  53.70975
#> 2 Sample 2 0.7    0.006  18.70133 3.574933 2.155759  18.68351  14.60380
#> 3 Sample 3 0.8    0.006 -10.28279 2.761595 1.739932 -10.27637 -13.73321
#>       upper
#> 1 66.796937
#> 2 23.005687
#> 3 -6.900347
```

All temperatures are in °C. `meanTemp` is the temperature from the mean
slope and intercept. `SE` is the standard error of the temperature,
based on `D47error`, the standard error of the mean D47. `SD` is the
uncertainty of the temperature of a single replicate, based on the SD of
the replicates (`D47error` × sqrt(`N`)). Both include the uncertainty in
the calibration. `SD` is `NA` when the dataset has neither `N` nor
`D47sd`.

`median`, `lower`, and `upper` are the median and the limits of the 95%
interval, based on `SE`. Temperature is a nonlinear function of D47, so
the interval is not symmetric around the median.

From here,
[`vignette("frequentist", package = "bayclumpr")`](https://tripati-lab.github.io/bayclumpr/articles/frequentist.md)
compares the four frequentist models.
[`vignette("bayesian", package = "bayclumpr")`](https://tripati-lab.github.io/bayclumpr/articles/bayesian.md)
shows the same workflow under a Bayesian framework.
