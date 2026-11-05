# Package index

## Prepare data

Convert temperatures, average replicates, and simulate datasets.

- [`temp2x()`](https://tripati-lab.github.io/bayclumpr/reference/temp2x.md)
  [`x2temp()`](https://tripati-lab.github.io/bayclumpr/reference/temp2x.md)
  : Convert temperatures to and from 10^6/T^2
- [`avg.replicates()`](https://tripati-lab.github.io/bayclumpr/reference/avg.replicates.md)
  : Average replicate measurements of D47
- [`cal.dataset()`](https://tripati-lab.github.io/bayclumpr/reference/cal.dataset.md)
  : Generate a synthetic dataset for clumped isotopes calibrations

## Calibrate (frequentist)

Bootstrap calibrations. Each returns a table of alpha, beta, and sigma.

- [`cal.ols()`](https://tripati-lab.github.io/bayclumpr/reference/cal.ols.md)
  : Fit OLS regression models on a given calibration dataset
- [`cal.wols()`](https://tripati-lab.github.io/bayclumpr/reference/cal.wols.md)
  : Fit weighted OLS regression models on a given calibration dataset
- [`cal.deming()`](https://tripati-lab.github.io/bayclumpr/reference/cal.deming.md)
  : Fit Deming regression models on a given calibration dataset
- [`cal.york()`](https://tripati-lab.github.io/bayclumpr/reference/cal.york.md)
  : Fit York regression models on a given calibration dataset

## Calibrate (Bayesian)

Calibrations in Stan, posterior draws, and model comparison.

- [`cal.bayesian()`](https://tripati-lab.github.io/bayclumpr/reference/cal.bayesian.md)
  [`print(`*`<bayclump_cal>`*`)`](https://tripati-lab.github.io/bayclumpr/reference/cal.bayesian.md)
  : Bayesian calibration of the clumped isotope paleothermometer
- [`cal.ion.bayesian()`](https://tripati-lab.github.io/bayclumpr/reference/cal.ion.bayesian.md)
  : Bayesian calibration with a second predictor
- [`cal.draws()`](https://tripati-lab.github.io/bayclumpr/reference/cal.draws.md)
  : Posterior draws from a Bayesian calibration
- [`cal.compare()`](https://tripati-lab.github.io/bayclumpr/reference/cal.compare.md)
  : Compare Bayesian calibration models

## Priors

Draw from the priors used in the Bayesian calibrations.

- [`cal.prior()`](https://tripati-lab.github.io/bayclumpr/reference/cal.prior.md)
  : Generate a dataset reflecting the priors used to run the analyses
- [`rec.prior()`](https://tripati-lab.github.io/bayclumpr/reference/rec.prior.md)
  : Generate a dataset reflecting the priors used to run the analyses

## Reconstruct

Reconstruct temperatures from a table of alpha, beta, and sigma.

- [`rec.clumped()`](https://tripati-lab.github.io/bayclumpr/reference/rec.clumped.md)
  : Temperature reconstructions from bootstrap calibrations
- [`rec.bayesian()`](https://tripati-lab.github.io/bayclumpr/reference/rec.bayesian.md)
  : Bayesian temperature reconstructions
- [`rec.ion.bayesian()`](https://tripati-lab.github.io/bayclumpr/reference/rec.ion.bayesian.md)
  : Bayesian temperature reconstructions with a second predictor

## Summarize and plot

Estimates with their standard errors, and calibration lines with their
intervals.

- [`cal.summary()`](https://tripati-lab.github.io/bayclumpr/reference/cal.summary.md)
  : Summarize a calibration
- [`cal.ci()`](https://tripati-lab.github.io/bayclumpr/reference/cal.ci.md)
  : Calibration line with its interval
