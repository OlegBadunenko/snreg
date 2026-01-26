# Compute Owen-like Function `TOwen1(h, a)`

`TOwen1` computes an Owen's \\T\\-function variant (or a related special
function) for vectors `h` and `a` based on the `tha` function in
<https://people.sc.fsu.edu/~jburkardt/c_src/owen/owen.html>. Non-finite
inputs in `h` or `a` yield `NA` at the corresponding positions.

## Usage

``` r
TOwen1(h, a, threads = 1)
```

## Arguments

- threads:

  integer. Number of threads to request from the C implementation (if
  supported). Default is `1`.

## Value

A numeric vector of length `length(h)` with the computed values.
Elements where either `h` or `a` is non-finite are `NA`. The returned
vector is given class `"snreg"` for downstream compatibility.

## Details

Owen's T Function Variant via C Backend

This is a thin R wrapper around a native routine with signature:


      void TOwen1(int *n, double *h, double *a, double *out, int *threads)

## See also

[`TOwen`](https://olegbadunenko.github.io/snreg/reference/TOwen.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  library(snreg)

  # Basic usage. Vectorized 'a':
  h <- c(-1, 0, 1, 2)
  a <- 0.3
  TOwen1(h, a)

  # Vectorized 'a' with non-finite entries:
  a2 <- c(0.2, NA, 1, Inf)
  TOwen1(h, a2)
} # }
```
