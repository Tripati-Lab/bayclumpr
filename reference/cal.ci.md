# This function is used to generate CI estimates at given intervals. It is currently used for plotting in BayClump.

This function is used to generate CI estimates at given intervals. It is
currently used for plotting in BayClump.

## Usage

``` r
cal.ci(data, from, to, length.out = 100)
```

## Arguments

- data:

  A `data.frame` with two columns named as beta and alpha. This should
  be the result of bootstrapping or the posterior distribution for a
  given calibration set.

- from:

  the lower limit in x.

- to:

  the upper limit in x.

- length.out:

  the number of breaks.
