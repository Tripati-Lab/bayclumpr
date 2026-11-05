# Generate a dataset reflecting the priors used to run the analyses

Same as
[`cal.prior()`](https://tripati-lab.github.io/bayclumpr/reference/cal.prior.md),
with two extra columns: `x`, values of 10^6/T^2 close to 11, and `xC`,
the same values in degrees C.

## Usage

``` r
rec.prior(prior, n = 1000)
```

## Arguments

- prior:

  The priors on alpha and beta. See
  [`cal.prior()`](https://tripati-lab.github.io/bayclumpr/reference/cal.prior.md).

- n:

  number of observations to simulate

## Value

A `data.frame` with `n` rows and the columns `alpha`, `beta`, `x`, and
`xC`.
