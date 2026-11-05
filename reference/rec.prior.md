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
#>         alpha       beta        x       xC
#> 1  0.36722189  1.5128812 11.00010 28.36003
#> 2  0.63816760  0.7162685 10.99994 28.36217
#> 3  0.16134519  0.4189627 10.99992 28.36238
#> 4 -0.01666434 -0.1537984 10.99984 28.36348
#> 5  0.92655081  1.6168918 10.99985 28.36334
#> 6  1.37722836  0.6352341 11.00001 28.36127
```
