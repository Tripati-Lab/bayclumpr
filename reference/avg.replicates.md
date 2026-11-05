# Average replicate measurements of D47

Averages the D47 replicates of each sample and returns the standard
error of each mean. The result can be used as a calibration or
reconstruction dataset.

## Usage

``` r
avg.replicates(data, pooledSD = NULL)
```

## Arguments

- data:

  A `data.frame` with one row per replicate and the columns `Sample` and
  `D47`.

- pooledSD:

  Optional SD used for every sample instead of the SD of its own
  replicates, such as the long-term SD of standards measured in the same
  lab.

## Value

A `data.frame` with one row per sample and the columns `Sample`, `D47`
(mean), `D47sd` (SD of the replicates), `N` (number of replicates), and
`D47error` (standard error of the mean), followed by the columns that
were kept.

## Details

`D47error` is the standard error of the mean, SD / sqrt(N). By default,
the SD comes from the replicates of each sample. Samples with a single
replicate then have no SE, so a warning is given. Other columns are kept
when they have a single value within every sample, such as `Material` or
`Temperature`.

## Examples

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
avg.replicates(reps, pooledSD = 0.008)
#>   Sample       D47       D47sd N    D47error Material
#> 1      A 0.6013333 0.003511885 3 0.004618802  calcite
#> 2      B 0.6986667 0.003511885 3 0.004618802  calcite
```
