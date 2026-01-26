# Summary for Skew-Normal Regression Models

Produces a summary object for objects of class `"snreg"`. The function
assigns the class `"summary.snreg"` to the fitted model object, enabling
a dedicated print method (`print.summary.snreg`) to display results in a
structured format.

## Usage

``` r
# S3 method for class 'snreg'
summary(obj, ...)
```

## Arguments

- obj:

  an object of class `"snreg"`, typically returned by
  [`snreg`](https://olegbadunenko.github.io/snreg/reference/snreg.md).

- ...:

  additional arguments (currently not used).

## Value

An object of class `"summary.snreg"`, identical to the input `obj`
except for its class attribute.

## Details

Summary Method for snreg Objects

This method expects a fitted `"snreg"` object.

`summary.snreg` does not modify the contents of the object; it only
updates the class attribute to `"summary.snreg"`. The corresponding
print method
([`print.summary.snreg`](https://olegbadunenko.github.io/snreg/reference/print.summary.snreg.md))
is responsible for formatting and displaying estimation details, such as
convergence criteria, log-likelihood, coefficient tables, and (if
present) heteroskedastic and skewness components.

## See also

[`snreg`](https://olegbadunenko.github.io/snreg/reference/snreg.md),
[`print.summary.snreg`](https://olegbadunenko.github.io/snreg/reference/print.summary.snreg.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  # m <- snreg(TC ~ Y1 + Y2, data = banks07)
  # s <- summary(m)
  # print(s)
} # }
```
