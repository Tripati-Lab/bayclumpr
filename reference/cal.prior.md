# Generate a dataset reflecting the priors used to run the analyses

Generate a dataset reflecting the priors used to run the analyses

## Usage

``` r
cal.prior(prior, n = 1000)
```

## Arguments

- prior:

  The priors on alpha and beta. Use `"Informative"` or `"Weak"`, a list
  such as `list(alpha = c(0.231, 0.065), beta = c(0.039, 0.004))`, or a
  result from
  [`cal.bayesian()`](https://tripati-lab.github.io/bayclumpr/reference/cal.bayesian.md)
  to use its priors. Flat priors are improper and can't be sampled.

- n:

  number of observations to simulate

## Value

A `data.frame` with `n` draws of `alpha` and `beta`.

## Examples

``` r
pr <- cal.prior("Informative", n = 1000)
quantile(pr$beta, c(0.025, 0.975))
#>       2.5%      97.5% 
#> 0.03084583 0.04666279 
cal.prior(list(alpha = c(0.25, 0.05), beta = c(0.038, 0.003)), n = 5)
#>       alpha       beta
#> 1 0.2881646 0.03730484
#> 2 0.2900645 0.03364383
#> 3 0.2714382 0.04091573
#> 4 0.2955264 0.03248460
#> 5 0.2320593 0.04042058
```
