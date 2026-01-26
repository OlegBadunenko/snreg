# Compute Owen's T Function `T(h, a)`

`TOwen1` computes an Owen's \\T\\-function variant (or a related special
function) for vectors `h` and `a` based on the `t` function in
<https://people.sc.fsu.edu/~jburkardt/c_src/owen/owen.html>. Non-finite
inputs (in `h` or `a`) produce `NA` at corresponding positions, while
finite pairs are computed in C in a vectorized fashion.

## Usage

``` r
TOwen(h, a, threads = 1)
```

## Arguments

- h:

  numeric vector of \\h\\ arguments.

- a:

  numeric vector of \\a\\ arguments. Must be either the same length as
  `h` or of length 1 (will be recycled by standard R rules).

- threads:

  integer. Number of threads to request from the C implementation (if
  supported). Default is `1`.

## Value

A numeric vector of length `length(h)` containing \\T(h_i, a_i)\\.
Elements where either `h_i` or `a_i` is not finite are `NA`. The
returned object is given class `"snreg"` for downstream compatibility
with your package’s print/summary helpers.

## Details

Owen's T Function via C Backend

Owen's \\T\\ function is commonly defined as \$\$T(h, a) \\=\\
\frac{1}{2\pi} \int\_{0}^{a} \frac{\exp\\\left(-\tfrac{1}{2}h^2
(1+t^2)\right)}{1+t^2} \\ dt,\$\$ for real \\h\\ and \\a\\.

The function accepts vector inputs and:

- Computes results only for entries where both `h` and `a` are finite.

- Returns `NA` where either `h` or `a` is non-finite.

- Optionally passes a `threads` hint to the C backend (ignored if not
  supported).

## See also

[`pnorm`](https://rdrr.io/r/stats/Normal.html),
[`dnorm`](https://rdrr.io/r/stats/Normal.html)

## Examples

``` r
if (FALSE) { # \dontrun{
  # Basic usage. Vectorized 'a'
  h <- c(-1, 0, 1, 2)
  a <- 0.5
  TOwen(h, a)

  # Vectorized 'a' with non-finite entries; non-finite entries yield NA
  a2 <- c(0.2, NA, 1, Inf)
  TOwen(h, a2)
} # }
```
