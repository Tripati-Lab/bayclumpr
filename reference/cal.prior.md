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
