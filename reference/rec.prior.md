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

## Examples

``` r
head(rec.prior("Weak", n = 100))
#>        alpha          beta        x       xC
#> 1 -1.3036063  1.3271487112 11.00001 28.36122
#> 2 -1.8196919 -1.0414222732 11.00004 28.36077
#> 3  1.0472887 -0.8799242990 10.99995 28.36202
#> 4 -0.1481982  0.7018517221 11.00003 28.36092
#> 5 -2.3567259  0.5516542666 10.99999 28.36144
#> 6  0.2189088  0.0005932573 10.99996 28.36185
```
