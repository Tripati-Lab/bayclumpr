<!-- badges: start -->
  [![R-CMD-check](https://github.com/Tripati-Lab/bayclumpr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Tripati-Lab/bayclumpr/actions/workflows/R-CMD-check.yaml)
[![CodeFactor](https://www.codefactor.io/repository/github/tripati-lab/bayclumpr/badge)](https://www.codefactor.io/repository/github/tripati-lab/bayclumpr)
  [![](https://img.shields.io/github/languages/code-size/tripati-lab/bayclumpr.svg)](https://github.com/tripati-lab/bayclumpr)
  [![](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html#maturing)  

<!-- badges: end -->


# The `bayclumpr` `R` Package

### Bayesian methods for clumped isotope paleothermometry


## What is `bayclumpr`?

To support the use of Bayesian models and the analytical framework developed in Román-Palacios et al. (2022) for clumped isotope calibration and for temperature reconstructions, and to facilitate comparisons of Bayesian and classical models, we present a self-contained `R` package and associated Shiny Dashboard application, `bayclumpr` and `BayClump`, respectively. `bayclumpr` (and `BayClump`) fits both frequentist and Bayesian linear regressions to calibration datasets and performs temperature reconstructions under both frameworks. You can find more details on how to use `bayclumpr` in this [website](https://tripati-lab.github.io/bayclumpr/).

## What is `BayClump`?

`BayClump` is an associated Shiny Dashboard application that is associated with the `bayclumpr` `R` package. All the functions implemented in `BayClump` are sourced from `bayclumpr`. We eveloped `BayClump` to allow users with less coding experience to be able to access standarized resources to calibrate and derive reconstructions using clumped isotope datasets. You can access BayClump directly from [here](https://bayclump.tripatilab.epss.ucla.edu/).    

## Installing `bayclumpr`

The latest version of `bayclumpr` can be installed directly from our GitHub repository using the following lines of code from `R`

```
library(devtools)
install_github("Tripati-Lab/bayclumpr")
```

## Contributing

Please see our [contributing guide](CONTRIBUTING).

## Contact

Please see the package [DESCRIPTION](DESCRIPTION) for package authors.

