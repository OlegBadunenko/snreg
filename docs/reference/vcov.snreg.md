# Extract the Variance–Covariance Matrix

`vcov.snreg` is the `vcov` S3 method for objects of class `"snreg"`. It
returns the model-based variance–covariance matrix stored in the fitted
object.

## Usage

``` r
# S3 method for class 'snreg'
vcov(obj, ...)
```

## Arguments

- obj:

  an object of class `"snreg"`, typically returned by
  [`snreg`](https://olegbadunenko.github.io/snreg/reference/snreg.md).

- ...:

  additional arguments (currently unused).

## Value

A numeric matrix containing the variance–covariance of the estimated
parameters.

## Details

Variance–Covariance Matrix for snreg Objects

This method expects a fitted `"snreg"` object.

This method simply returns the `vcov` component stored in `obj`. If your
estimator did not compute standard errors (e.g., because estimation
hasn’t been run yet in a scaffold), this field may be `NULL`, and the
method will error accordingly.

## See also

[`snreg`](https://olegbadunenko.github.io/snreg/reference/snreg.md),
[`summary.snreg`](https://olegbadunenko.github.io/snreg/reference/summary.snreg.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  library(snsf)

  data("banks07")
  head(banks07)
  # V <- vcov(m)
  # diag(V)
} # }
```
