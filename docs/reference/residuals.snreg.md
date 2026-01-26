# Residuals for snreg Objects

`residuals.snreg` is the S3 method for extracting residuals from a
fitted `snreg` model. Residuals may be returned either for the full data
or only for the estimation sample.

## Usage

``` r
# S3 method for class 'snreg'
residuals(obj, esample = TRUE, ...)
```

## Arguments

- obj:

  an object of class `"snreg"`, typically produced by
  [`snreg`](https://olegbadunenko.github.io/snreg/reference/snreg.md).

- esample:

  logical. If `TRUE` (default), residuals are returned only for
  observations used in estimation (others are `NA`). If `FALSE`, the raw
  vector of residuals (`obj$resid`) is returned.

- ...:

  additional arguments (currently unused).

## Value

A numeric vector of residuals. If `esample = TRUE`, the vector matches
the length of the original data and contains `NA` for non-estimation
observations. If `esample = FALSE`, only the computed residuals are
returned.

## Details

Extract Residuals from an snreg Model

This method simply accesses the `obj$resid` component of a fitted
`"snreg"` object. An informative error is produced if residuals are not
available.

## See also

[`snreg`](https://olegbadunenko.github.io/snreg/reference/snreg.md),
`fitted.snreg`,
[`coef.snreg`](https://olegbadunenko.github.io/snreg/reference/coef.snreg.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  m <- snreg(y ~ x1 + x2, data = df)

  # Residuals for estimation sample only
  residuals(m)

  # Residuals for all observations
  residuals(m, esample = FALSE)
} # }
```
