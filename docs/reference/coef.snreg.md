# Extract Model Coefficients

`coef.snreg` is the S3 method for extracting the estimated regression
coefficients from an object of class `"snreg"`.

## Usage

``` r
# S3 method for class 'snreg'
coef(obj, ...)
```

## Arguments

- obj:

  an object of class `"snreg"`, typically returned by
  [`snreg`](https://olegbadunenko.github.io/snreg/reference/snreg.md).

- ...:

  additional arguments (currently unused).

## Value

A numeric vector containing the model coefficients.

## Details

Coefficients from an snreg Model

This method simply returns the `coef` component stored inside the fitted
`"snreg"` object. If the object does not contain coefficient estimates
(e.g., if estimation was not completed in a scaffold), an informative
error is raised.

## Examples

``` r
if (FALSE) { # \dontrun{
  m <- snreg(y ~ x1 + x2, data = df)
  coef(m)
} # }
```
