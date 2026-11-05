# Package index

## All functions

- [`cal.bayesian()`](https://tripati-lab.github.io/bayclumpr/reference/cal.bayesian.md)
  :

  Bayesian regressions to calibrate the clumped isotopes
  paleothermometer using `stan`.

- [`cal.ci()`](https://tripati-lab.github.io/bayclumpr/reference/cal.ci.md)
  : This function is used to generate CI estimates at given intervals.
  It is currently used for plotting in BayClump.

- [`cal.dataset()`](https://tripati-lab.github.io/bayclumpr/reference/cal.dataset.md)
  : Generate a synthetic dataset for clumped isotopes calibrations

- [`cal.deming()`](https://tripati-lab.github.io/bayclumpr/reference/cal.deming.md)
  : Fit Deming regression models on a given calibration dataset

- [`cal.ion.bayesian()`](https://tripati-lab.github.io/bayclumpr/reference/cal.ion.bayesian.md)
  :

  Bayesian regressions to calibrate the clumped isotopes
  paleothermometer using `stan`.

- [`cal.ols()`](https://tripati-lab.github.io/bayclumpr/reference/cal.ols.md)
  : Fit OLS regression models on a given calibration dataset

- [`cal.prior()`](https://tripati-lab.github.io/bayclumpr/reference/cal.prior.md)
  : Generate a dataset reflecting the priors used to run the analyses

- [`cal.wols()`](https://tripati-lab.github.io/bayclumpr/reference/cal.wols.md)
  : Fit weighted OLS regression models on a given calibration dataset

- [`cal.york()`](https://tripati-lab.github.io/bayclumpr/reference/cal.york.md)
  : Fit York regression models on a given calibration dataset

- [`rec.bayesian()`](https://tripati-lab.github.io/bayclumpr/reference/rec.bayesian.md)
  : This function generate temperature predictions (in 10^6/T2) based on
  a calibration dataset and target D47. Note that this approach
  additionally accounts for measured error in the target D47. This
  approach is congruent with the one used in McClelland et al. (2022).

- [`rec.clumped()`](https://tripati-lab.github.io/bayclumpr/reference/rec.clumped.md)
  : This function performs temp reconstruction (10^6/T^2 with T in K)
  for multiple replicates of the same target.

- [`rec.ion.bayesian()`](https://tripati-lab.github.io/bayclumpr/reference/rec.ion.bayesian.md)
  : This function generates temperature predictions (in 10^6/T2) based
  on a calibration dataset and target D47. It includes measured error in
  the target D47 and ion as an additional predictor.

- [`rec.prior()`](https://tripati-lab.github.io/bayclumpr/reference/rec.prior.md)
  : Generate a dataset reflecting the priors used to run the analyses
