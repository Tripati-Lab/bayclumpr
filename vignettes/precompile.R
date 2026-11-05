# Vignettes that fit Stan models are precompiled so the package builds quickly
# Run this script from the package root after changing bayesian.Rmd.orig or covariates.Rmd.orig

library(knitr)
knit("vignettes/bayesian.Rmd.orig", "vignettes/bayesian.Rmd")
knit("vignettes/covariates.Rmd.orig", "vignettes/covariates.Rmd")

library(devtools)
build_vignettes()
