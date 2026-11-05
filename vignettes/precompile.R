# Vignettes that fit Stan models are precompiled so the package builds quickly
# Run this script from the package root after changing any .Rmd.orig file
# Knitting from the vignettes folder keeps the figures next to the vignettes

library(knitr)
old <- setwd("vignettes")
knit("bayesian.Rmd.orig", "bayesian.Rmd")
knit("covariates.Rmd.orig", "covariates.Rmd")
knit("comparison.Rmd.orig", "comparison.Rmd")
setwd(old)

library(devtools)
build_vignettes()
