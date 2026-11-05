# This function generate temperature predictions (in 10^6/T2) based on a calibration dataset and target D47. Note that this approach additionally accounts for measured error in the target D47. This approach is congruent with the one used in McClelland et al. (2022).

This function generate temperature predictions (in 10^6/T2) based on a
calibration dataset and target D47. Note that this approach additionally
accounts for measured error in the target D47. This approach is
congruent with the one used in McClelland et al. (2022).

## Usage

``` r
rec.bayesian(
  calModel,
  recData,
  iter = 1000,
  mixed = FALSE,
  postcalsamples = NULL,
  MC = TRUE
)
```

## Arguments

- calModel:

  The stan model to be analyzed.

- recData:

  The reconstruction dataset.

- iter:

  Number of replicates to retain.

- mixed:

  whether the model `calModel` is mixed or not.

- postcalsamples:

  Number of posterior samples to analyze from the calibration step.

- MC:

  Multicore (TRUE/FALSE)
