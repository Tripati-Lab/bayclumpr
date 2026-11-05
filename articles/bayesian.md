# 3. Bayesian calibrations and reconstructions

This vignette covers Bayesian calibrations in `bayclumpr`. We fit the
models, choose priors, check convergence, compare models, and
reconstruct temperatures. If you’re new to the package, start with
[`vignette("bayclumpr", package = "bayclumpr")`](https://tripati-lab.github.io/bayclumpr/articles/bayclumpr.md).

``` r

library(bayclumpr)
```

## The four models

[`cal.bayesian()`](https://tripati-lab.github.io/bayclumpr/reference/cal.bayesian.md)
fits up to four models. They cover every combination of two choices:
whether measurement errors are used, and whether intercepts and slopes
vary by `Material`. We choose them with the `models` argument.

| `models` | Name in the result | Measurement errors | Varies by `Material` |
|:---|:---|:--:|:--:|
| `"errors"` | `BLM1_fit` | Yes | No |
| `"noerrors"` | `BLM1_fit_NoErrors` | No | No |
| `"mixed_errors"` | `BLM3_fit_Errors` | Yes | Yes |
| `"mixed_noerrors"` | `BLM3_fit` | No | Yes |

Models with measurement errors use `TempError` and `D47error` for each
observation. All models estimate `sigma`, the scatter around the
calibration line. With measurement errors, `sigma` is the scatter left
after removing the measurement error in `D47`. Without them, it includes
the measurement error.

Mixed models let intercepts and slopes vary by `Material` around shared
means. They need at least two materials and are skipped otherwise.

The models run in Stan. If they weren’t compiled when `bayclumpr` was
installed, the first fit of each model in a session compiles it. This
can take a minute or more. Later fits in the same session start right
away.

## Fitting a model

We simulate 50 observations under the intermediate error scenario. The
true values are `alpha` = 0.268 and `beta` = 0.0369.

``` r

set.seed(3)
ds <- cal.dataset(error = "S2", nobs = 50)
```

We start with the errors model:

``` r

fit <- cal.bayesian(calibrationData = ds, models = "errors")
fit
#> Bayesian clumped isotope calibration
#> Priors: alpha ~ N(0.231, 0.065), beta ~ N(0.039, 0.004)
#> 
#> BLM1_fit
#>  parameter estimate       SE  median   lower   upper
#>      alpha  0.25881 0.020695 0.25869 0.21710 0.29843
#>       beta  0.03750 0.001717 0.03749 0.03419 0.04095
#>      sigma  0.02463 0.002906 0.02439 0.01948 0.03066
#> 
#> Convergence
#>     model max_Rhat min_ESS divergences
#>  BLM1_fit   0.9997    3997           0
```

The first table is the same summary given by
[`cal.summary()`](https://tripati-lab.github.io/bayclumpr/reference/cal.summary.md).
`estimate` is the posterior mean of each parameter and `SE` is its
standard error, which is the SD of the posterior draws. `sigma` is
itself an SD: the scatter of individual samples around the calibration
line. The second table checks convergence. We look for `max_Rhat` at or
below 1.01, `min_ESS` of at least 100 per chain, and no divergences.
[`cal.bayesian()`](https://tripati-lab.github.io/bayclumpr/reference/cal.bayesian.md)
gives a warning when a model fails any of these checks. When that
happens, try a larger `numSavedSteps` and check the data for outliers or
wrong units.

## Choosing priors

The `priors` argument sets the priors on `alpha` and `beta`. There are
three presets:

| `priors`        | `alpha`              | `beta`               |
|:----------------|:---------------------|:---------------------|
| `"Informative"` | Normal(0.231, 0.065) | Normal(0.039, 0.004) |
| `"Weak"`        | Normal(0.231, 1)     | Normal(0.039, 1)     |
| `"Flat"`        | Flat                 | Flat                 |

The informative values are the ones used in earlier versions of this
package. Calibrations differ among reference frames and acid
corrections, so check that these values suit your data. You can also set
your own priors as a list with the mean and SD of each normal prior:

``` r

myPriors <- list(alpha = c(0.25, 0.05), beta = c(0.038, 0.003))
```

[`cal.prior()`](https://tripati-lab.github.io/bayclumpr/reference/cal.prior.md)
draws from a set of priors. This shows the range of slopes that a prior
allows before seeing the data:

``` r

quantile(cal.prior("Informative", n = 10000)$beta, c(0.025, 0.975))
#>       2.5%      97.5% 
#> 0.03105497 0.04668649
quantile(cal.prior(myPriors, n = 10000)$beta, c(0.025, 0.975))
#>       2.5%      97.5% 
#> 0.03203628 0.04367196
```

Let’s fit the same model under weak priors and compare the slopes:

``` r

fitWeak <- cal.bayesian(calibrationData = ds, models = "errors", priors = "Weak")

rbind(
  Informative = quantile(cal.draws(fit)$beta, c(0.025, 0.5, 0.975)),
  Weak = quantile(cal.draws(fitWeak)$beta, c(0.025, 0.5, 0.975))
)
#>                   2.5%        50%      97.5%
#> Informative 0.03419203 0.03748993 0.04094567
#> Weak        0.03260000 0.03677161 0.04075844
```

## Comparing models

We can fit several models at once and compare them with
[`cal.compare()`](https://tripati-lab.github.io/bayclumpr/reference/cal.compare.md).
It uses approximate leave-one-out cross-validation.

``` r

fits <- cal.bayesian(calibrationData = ds, models = c("errors", "noerrors"))
cal.compare(fits)
#>              model elpd_diff se_diff p_worse diag_diff diag_elpd
#>  BLM1_fit_NoErrors       0.0     0.0      NA                    
#>           BLM1_fit      -0.1     0.1    0.89   N < 100
#> 
#> Diagnostic flags present.
#> See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
#> or https://mc-stan.org/loo/reference/loo-glossary.html.
```

The best model is in the first row. When `elpd_diff` is small relative
to `se_diff`, the data don’t clearly favor one model.

## Mixed models

To show the mixed models, we simulate a second material with D47 values
that are 0.02 higher:

``` r

ds1 <- cal.dataset(error = "S2", nobs = 50)
ds2 <- cal.dataset(error = "S2", nobs = 50)
ds2$D47 <- ds2$D47 + 0.02
ds2$Material <- 2
dsMixed <- rbind(ds1, ds2)

fitMixed <- cal.bayesian(calibrationData = dsMixed, models = c("mixed_errors", "mixed_noerrors"))
#> Warning: There were 24 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: There were 23 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> Warning in .finish_fit(fits, prior, chains): Check convergence for: BLM3_fit_Errors,
#> BLM3_fit. Print the result for details.
fitMixed
#> Bayesian clumped isotope calibration
#> Priors: alpha ~ N(0.231, 0.065), beta ~ N(0.039, 0.004)
#> 
#> BLM3_fit_Errors
#>  Material parameter estimate       SE  median   lower   upper
#>         1     alpha  0.26988 0.016214 0.26962 0.23926 0.30229
#>         1      beta  0.03747 0.001399 0.03748 0.03463 0.04016
#>         1     sigma  0.02435 0.001956 0.02421 0.02084 0.02837
#>         2     alpha  0.26863 0.021463 0.26867 0.22481 0.31017
#>         2      beta  0.03833 0.001765 0.03833 0.03483 0.04192
#>         2     sigma  0.02435 0.001956 0.02421 0.02084 0.02837
#> 
#> BLM3_fit
#>  Material parameter estimate       SE  median   lower   upper
#>         1     alpha  0.27035 0.016371 0.26990 0.23813 0.30256
#>         1      beta  0.03744 0.001407 0.03745 0.03462 0.04016
#>         1     sigma  0.02570 0.001887 0.02558 0.02224 0.02966
#>         2     alpha  0.26913 0.021344 0.26930 0.22539 0.31011
#>         2      beta  0.03828 0.001749 0.03827 0.03488 0.04184
#>         2     sigma  0.02570 0.001887 0.02558 0.02224 0.02966
#> 
#> Convergence
#>            model max_Rhat min_ESS divergences
#>  BLM3_fit_Errors    1.002    2692          24
#>         BLM3_fit    1.003    1739          23
```

The intercepts and slopes of the two materials are pulled toward shared
means (partial pooling). The SDs among materials have half-normal priors
with scales of 0.1 for the intercept at the mean calibration temperature
and 0.01 for the slope. With fewer than five materials, these SDs mostly
reflect their priors.

## Reconstructing temperatures

[`cal.draws()`](https://tripati-lab.github.io/bayclumpr/reference/cal.draws.md)
returns the posterior draws of a model as a table. It has the same
columns as the frequentist tables in
[`vignette("frequentist", package = "bayclumpr")`](https://tripati-lab.github.io/bayclumpr/articles/frequentist.md).

``` r

draws <- cal.draws(fit)
head(draws)
#>       alpha       beta      sigma
#> 1 0.2913857 0.03403879 0.02635096
#> 2 0.2286358 0.03961537 0.02440160
#> 3 0.2765121 0.03592458 0.02101763
#> 4 0.2381617 0.03879429 0.02470852
#> 5 0.2492998 0.03765331 0.01982520
#> 6 0.2539153 0.03745754 0.02471296
```

We create a reconstruction dataset and pass it to
[`rec.bayesian()`](https://tripati-lab.github.io/bayclumpr/reference/rec.bayesian.md)
with the draws:

``` r

recData <- data.frame(
  Sample = paste("Sample", 1:6),
  D47 = rep(c(0.6, 0.7, 0.8), 2),
  D47error = rep(c(0.005, 0.02), each = 3),
  N = 3
)

rec.bayesian(calModel = draws, recData = recData)
#>     Sample D47 D47error  meanTemp        SD        SE     median      lower     upper
#> 1 Sample 1 0.6    0.005 59.169525 13.492686 12.919351  58.141071  36.176558 86.901763
#> 2 Sample 2 0.7    0.005 18.679092  8.878960  8.507950  18.238345   2.936739 36.380216
#> 3 Sample 3 0.8    0.005 -9.781490  6.594663  6.322624  -9.933485 -21.692635  3.210007
#> 4 Sample 4 0.6    0.020 59.489531 21.596755 16.403661  58.415102  30.709281 95.152623
#> 5 Sample 5 0.7    0.020 18.852540 14.311770 10.866998  18.135342  -0.351929 42.237888
#> 6 Sample 6 0.8    0.020 -9.734762 10.696587  7.841322 -10.021680 -24.568806  5.917603
```

`meanTemp` is the posterior mean temperature in °C. `SE` is the standard
error of the temperature, based on `D47error`, the standard error of the
mean D47. `SD` is the uncertainty of the temperature of a single
replicate, based on `D47error` × sqrt(`N`). `SD` is `NA` when the
dataset has neither `N` nor `D47sd`. `median`, `lower`, and `upper`
describe the 95% posterior interval, based on `SE`. Both `SD` and `SE`
include the calibration, the scatter around the calibration line
(`sigma`), and the error in `D47`.
[`rec.clumped()`](https://tripati-lab.github.io/bayclumpr/reference/rec.clumped.md)
gives the same kind of uncertainty with `prediction = TRUE`.

[`rec.bayesian()`](https://tripati-lab.github.io/bayclumpr/reference/rec.bayesian.md)
doesn’t run another MCMC. For each posterior draw, the temperature has a
known distribution, so it’s sampled directly. This takes less than a
second.

We can also pass the result of
[`cal.bayesian()`](https://tripati-lab.github.io/bayclumpr/reference/cal.bayesian.md)
and name the model:

``` r

rec.bayesian(calModel = fits, recData = recData, model = "noerrors")
#>     Sample D47 D47error  meanTemp        SD        SE     median       lower     upper
#> 1 Sample 1 0.6    0.005 59.582528 13.850920 13.356362  58.760144  35.6641789 87.739793
#> 2 Sample 2 0.7    0.005 18.738866  9.209379  8.807402  18.453100   2.3154928 36.871627
#> 3 Sample 3 0.8    0.005 -9.680912  6.965055  6.621272  -9.908738 -22.0689063  3.782275
#> 4 Sample 4 0.6    0.020 59.743116 22.361167 16.759462  58.531780  30.6149554 95.806125
#> 5 Sample 5 0.7    0.020 19.093028 14.569364 10.940762  18.535532  -0.3365017 42.545939
#> 6 Sample 6 0.8    0.020 -9.661798 10.738637  8.158392 -10.011537 -24.7372205  6.941973
```

### Reconstructions with mixed models

For mixed models, the reconstruction dataset needs a `Material` column
with the same labels used in the calibration:

``` r

recMixed <- data.frame(
  Sample = c("Sample 1", "Sample 2"),
  D47 = c(0.62, 0.62),
  D47error = 0.01,
  Material = c(1, 2)
)

rec.bayesian(calModel = fitMixed, recData = recMixed, model = "mixed_errors")
#>     Sample  D47 D47error meanTemp SD       SE   median    lower    upper
#> 1 Sample 1 0.62     0.01 54.89756 NA 12.72386 54.08777 32.33342 82.54673
#> 2 Sample 2 0.62     0.01 58.08932 NA 13.09937 57.30202 34.82513 86.17134
```

Mixed models can also reconstruct temperatures for a material that is
not in the calibration. With `newMaterial = TRUE`, the intercept and
slope are drawn from the shared means and the spread among materials.
The intervals are wider, since we don’t know where the new material
falls among the others.

``` r

recNew <- data.frame(Sample = "Sample 3", D47 = 0.62, D47error = 0.01)
rec.bayesian(calModel = fitMixed, recData = recNew, model = "mixed_errors", newMaterial = TRUE)
#> Warning in .rec_exact(draws, recData, ndraws, level = level): 27 draws had 10^6/T^2 at or
#> below zero and were removed.
#>     Sample  D47 D47error meanTemp SD       SE   median    lower    upper
#> 1 Sample 3 0.62     0.01 60.44808 NA 62.03957 54.12968 9.942506 136.8542
```

## Tips for faster fits

- Fit only the models you need with the `models` argument.
- Run chains in parallel with `cores`, for example `cores = 4`.
- Use fewer draws while you explore, for example `numSavedSteps = 1000`,
  and more for final results.
