#Vignettes that depend on internet access have been pre-compiled:

library(knitr)
knit("vignettes/bayclumpr.Rmd.orig", "vignettes/bayclumpr.Rmd")

library(devtools)
build_vignettes()
