# Print Summary of snreg Results

Prints the contents of a `"summary.snreg"` object in a structured
format. The method reports convergence status (based on gradient-Hessian
scaling), log-likelihood, estimation results, and—when present—summaries
for technical/cost efficiencies and marginal effects.

## Usage

``` r
# S3 method for class 'summary.snreg'
print(obj, digits = NULL, ...)
```

## Arguments

- obj:

  an object of class `"summary.snreg"` (produced by
  [`summary.snreg`](https://olegbadunenko.github.io/snreg/reference/summary.snreg.md)).

- digits:

  integer indicating the number of digits to print; default `NULL`
  (internally set to 4).

- ...:

  additional arguments (currently unused).

## Value

The input `obj` is returned (invisibly) after printing.

## Details

Print Method for Summary of snreg Objects

This method expects a fitted `"snreg"` object.

## See also

[`summary.snreg`](https://olegbadunenko.github.io/snreg/reference/summary.snreg.md)
