# Convert temperatures to and from 10^6/T^2

`temp2x()` converts temperatures in degrees C or kelvin to 10^6/T^2 (T
in K), the scale used by all calibrations in `bayclumpr`. `x2temp()`
converts them back.

## Usage

``` r
temp2x(temp, tempError = NULL, unit = c("C", "K"))

x2temp(x, xError = NULL, unit = c("C", "K"))
```

## Arguments

- temp:

  Temperatures in degrees C or kelvin.

- tempError:

  Optional 1 SD of `temp`, in the same units.

- unit:

  Unit of `temp`, or the unit returned by `x2temp()`: `"C"` or `"K"`.

- x:

  Temperatures as 10^6/T^2 (T in K).

- xError:

  Optional 1 SD of `x`.

## Value

`temp2x()` returns a `data.frame` with `Temperature` and, when
`tempError` is given, `TempError`. These names match the columns needed
by the calibration functions. `x2temp()` returns a `data.frame` with
`Temp` and, when `xError` is given, `TempSD`.

## Details

Errors are propagated to first order. The SD of 10^6/T^2 is 2 x 10^6 x
SD(T) / T^3, with T in kelvin. This approximation works well when the
error is small relative to T, which is the case for calibration data.

## Examples

``` r
temp2x(c(5, 25, 50), tempError = c(0.5, 0.5, 1))
#>   Temperature  TempError
#> 1   12.925337 0.04646894
#> 2   11.249426 0.03773076
#> 3    9.576166 0.05926762
x2temp(c(12.9, 11.3, 9.6), xError = 0.02)
#>        Temp    TempSD
#> 1  5.273023 0.2158318
#> 2 24.332059 0.2632585
#> 3 49.598612 0.3361965
```
