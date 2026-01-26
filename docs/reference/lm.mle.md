# Linear Regression via MLE

`lm.mle` fits a linear regression model by maximum likelihood, allowing
for optional multiplicative heteroskedasticity in the disturbance
variance via a log-linear specification provided through `ln.var.v`.

## Usage

``` r
lm.mle(
  formula,
  data,
  subset,
  ln.var.v = NULL,
  technique = c("bfgs"),
  lmtol = 1e-05,
  reltol = 1e-12,
  maxit = 199,
  optim.report = 1,
  optim.trace = 10,
  print.level = 3,
  digits = 4,
  only.data = FALSE,
  ...
)
```

## Arguments

- formula:

  an object of class `formula` specifying the regression: typically
  `y ~ x1 + ...`, where `y` is the dependent variable and the `x`'s are
  regressors.

- data:

  an optional `data.frame` containing the variables referenced in
  `formula`. If not found in `data`, variables are taken from
  `environment(formula)`.

- subset:

  an optional logical or numeric vector specifying the subset of
  observations to be used in estimation.

- ln.var.v:

  optional one-sided formula; e.g. `ln.var.v ~ z1 + z2`. When provided,
  the error variance is modeled as \\\log(\sigma_i^2) = w_i^\top
  \gamma_v\\. If `NULL`, the variance is homoskedastic.

- technique:

  character vector specifying the preferred optimization routine(s) in
  order of preference. Recognized keywords (for future implementation)
  include `"bfgs"` `"bhhh"`, `"nm"` (Nelder–Mead), `"bfgs"`, and `"cg"`.
  Default is `"bfgs"`. This scaffold records but does not execute the
  chosen routine.

- lmtol:

  numeric. Convergence tolerance based on scaled gradient (when
  applicable). Default `1e-5`.

- reltol:

  numeric. Relative convergence tolerance for likelihood maximization.
  Default `1e-12`.

- maxit:

  integer. Maximum number of iterations for the optimizer. Default
  `199`.

- optim.report:

  integer. Verbosity level for reporting progress (if implemented).
  Default `1`.

- optim.trace:

  integer. Trace level for optimization (if implemented). Default `1`.

- print.level:

  integer. Printing level for summaries. Default `3`.

- digits:

  integer. Number of digits for printing. Default `4`.

- only.data:

  logical. If `TRUE`, returns only constructed data/matrices without
  estimation. Default `FALSE`.

- ...:

  additional arguments reserved for future methods (e.g., bounds,
  penalties).

## Value

A list of class `"snreg"` containing (and extending) the fields returned
by [`optim`](https://rdrr.io/r/stats/optim.html):

- `par` — numeric vector of the MLE parameter estimates.

- `value` — numeric scalar: maximized log-likelihood value.

- `ll` — numeric scalar: maximized log-likelihood value.

- `counts`, `convergence`, `message` — standard `optim` outputs.

- `hessian` — the observed Hessian at the solution (as returned by
  `optim(hessian=TRUE)`).

- `coef` — named numeric vector equal to `par` (estimates).

- `vcov` — variance–covariance matrix, computed as `solve(-hessian)`.

- `sds` — standard errors, `sqrt(diag(vcov))`.

- `ctab` — coefficient table with columns: `Estimate`, `Std.Err`,
  `Z value`, `Pr(>z)`.

- `esample` — logical vector indicating the observations used in
  estimation.

- `n` — scalar, number of observations in the estimation sample.

The object inherits the default `optim` components and is assigned class
`"snreg"`.

## Details

Linear Model by Maximum Likelihood (with optional heteroskedasticity)

The model is \$\$y_i = x_i^\top \beta + \varepsilon_i,\quad
\varepsilon_i \sim \mathcal{N}(0, \sigma_i^2).\$\$ When `ln.var.v` is
supplied, the variance follows \$\$\log(\sigma_i^2) = w_i^\top
\gamma_v,\$\$ otherwise \\\sigma_i^2 = \sigma^2\\ is constant
(homoskedastic).

This function:

- Builds the model frame and `X`, `y`.

- Builds `Zv` for the log-variance index when `ln.var.v` is provided.

- Returns a structured object with placeholders for `coef`, `vcov`,
  `loglik`.

Insert your MLE engine to estimate \\\beta\\, and (optionally)
\\\sigma^2\\ or \\\gamma_v\\; compute standard errors via AIM/OPG as
required by `vcetype`.

## Examples

``` r
if (FALSE) { # \dontrun{

library(snreg)

data("banks07")
head(banks07)

# Translog cost function specification
spe.tl <- log(TC) ~ (log(Y1) + log(Y2) + log(W1) + log(W2))^2 +
  I(0.5 * log(Y1)^2) + I(0.5 * log(Y2)^2) +
  I(0.5 * log(W1)^2) + I(0.5 * log(W2)^2)

# -------------------------------------------------------------
# Specification 1: homoskedastic noise (ln.var.v = NULL)
# -------------------------------------------------------------
formSV <- NULL

m1 <- lm.mle(
  formula   = spe.tl,
  data      = banks07,
  ln.var.v  = formSV
)

coef(m1)


# -------------------------------------------------------------
# Specification 2: heteroskedastic noise (variance depends on TA)
# -------------------------------------------------------------
formSV <- ~ log(TA)

m2 <- lm.mle(
  formula   = spe.tl,
  data      = banks07,
  ln.var.v  = formSV
)

coef(m2)

} # }
```
