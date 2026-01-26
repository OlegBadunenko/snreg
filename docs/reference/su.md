# Summary Utility (`su`)

Computes a compact table of summary statistics for each variable in a
vector, matrix, or data frame. The following metrics are returned per
variable: number of observations (`Obs`), missing values (`NAs`), mean,
standard deviation (`StDev`), interquartile range (`IQR`), minimum
(`Min`), user-specified quantiles (`probs`), and maximum (`Max`).

## Usage

``` r
su(
  x,
  mat.var.in.col = TRUE,
  digits = 4,
  probs = c(0.1, 0.25, 0.5, 0.75, 0.9),
  print = FALSE
)
```

## Arguments

- x:

  a numeric vector, matrix, or data frame. For matrices, variables are
  assumed to be in columns; set `mat.var.in.col = FALSE` to treat rows
  as variables.

- mat.var.in.col:

  logical. If `TRUE` (default), a matrix is interpreted as variables in
  columns. If `FALSE`, the matrix is transposed so that rows are treated
  as variables.

- digits:

  integer. Number of digits to use when printing (only affects printed
  output when `print = TRUE`). Default is `4`.

- probs:

  numeric vector of probabilities in \\\[0, 1\]\\ for which quantiles
  are computed. Default is `c(0.1, 0.25, 0.5, 0.75, 0.9)`.

- print:

  logical. If `TRUE`, prints the transposed summary table using the
  specified number of digits. Default is `FALSE`.

## Value

A matrix (coercible to `data.frame`) where each row corresponds to a
variable and columns contain the summary statistics: `Obs`, `NAs`,
`Mean`, `StDev`, `IQR`, `Min`, the requested `probs` quantiles (named),
and `Max`. The returned object is given class `"snreg"` for
compatibility with package-specific print/summarization methods.

## Details

Compact Summary Statistics for Vectors, Matrices, and Data Frames

Input handling:

- If `x` is a matrix with a single row or column, it is treated like a
  vector. Column or row names are used (if available). Otherwise, a
  default name is created.

- If `x` is a matrix with multiple variables, variables are taken as
  columns. Use `mat.var.in.col = FALSE` to transpose and treat rows as
  variables.

- If `x` is a vector, its deparsed symbol name is used as the variable
  name.

- If `x` is a data frame, each column is summarized.

Missing values are excluded in all summary computations.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Vector
  set.seed(1)
  v <- rnorm(100)
  su(v, print = TRUE)

  # Matrix: variables in columns
  M <- cbind(x = rnorm(50), y = runif(50))
  su(M)

  # Matrix: variables in rows
  Mr <- rbind(x = rnorm(50), y = runif(50))
  su(Mr, mat.var.in.col = FALSE)

  # Data frame
  DF <- data.frame(a = rnorm(30), b = rexp(30), c = rbinom(30, 1, 0.3))
  out <- su(DF)
  head(out)
} # }
```
