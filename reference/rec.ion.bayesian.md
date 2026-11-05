# This function generates temperature predictions (in 10^6/T2) based on a calibration dataset and target D47. It includes measured error in the target D47 and ion as an additional predictor.

This function generates temperature predictions (in 10^6/T2) based on a
calibration dataset and target D47. It includes measured error in the
target D47 and ion as an additional predictor.

## Usage

``` r
rec.ion.bayesian(
  calModel,
  recData,
  iter = 1000,
  postcalsamples = NULL,
  MC = TRUE,
  useIonError = FALSE
)
```

## Arguments

- calModel:

  The fitted Stan model from calibration.

- recData:

  The reconstruction dataset.

- iter:

  Number of sampling iterations.

- postcalsamples:

  Number of posterior calibration samples to use.

- MC:

  Enable multicore (default TRUE).

- useIonError:

  Use Ion error or not.
