# Stochastic Frontier Model with a Skew-Normally Distributed Error Term

`snsf` performs maximum likelihood estimation of the parameters and
technical or cost efficiencies in a Stochastic Frontier Model with a
skew-normally distributed error term.

## Usage

``` r
snsf(
  formula,
  data,
  subset,
  distribution = "e",
  prod = TRUE,
  start.val = NULL,
  init.sk = NULL,
  ln.var.u = NULL,
  ln.var.v = NULL,
  skew.v = NULL,
  mean.u = NULL,
  technique = c("nr"),
  vcetype = c("aim"),
  optim.report = 1,
  optim.trace = 1,
  reltol = 1e-12,
  lmtol = 1e-05,
  maxit = 199,
  print.level = 3,
  report = 1,
  trace = 1,
  threads = 1,
  only.data = FALSE,
  digits = 4,
  ...
)
```

## Arguments

- formula:

  an object of class `formula` specifying the frontier: a typical model
  is `y ~ x1 + ...`, where `y` is the log of output (or total cost), and
  `x`'s are inputs (or outputs and input prices, in logs). See
  **Details**.

- data:

  an optional `data.frame` containing the variables in `formula`. If not
  found in `data`, variables are taken from `environment(formula)`.

- subset:

  an optional logical or numeric vector specifying a subset of
  observations for which the model is estimated and efficiencies are
  computed.

- distribution:

  character scalar specifying the distribution of the inefficiency term:
  default `"e"` (exponential). `"h"` (half-normal) and `"t"` (truncated
  normal) to be implemented.

- prod:

  logical. If `TRUE`, estimates correspond to a stochastic *production*
  frontier and technical efficiencies are returned; if `FALSE`,
  estimates correspond to a stochastic *cost* frontier and cost
  efficiencies are returned. Default is `TRUE`.

- start.val:

  optional numeric vector of starting values for the optimizer.

- init.sk:

  numeric. Initial value for the skewness parameter of the noise
  component; default is `0.5`.

- ln.var.u:

  optional one-sided formula; e.g. `ln.var.u = ~ z3 + z4`. Specifies
  exogenous variables entering the (log) variance of the inefficiency
  component. If `NULL`, the inefficiency variance is homoskedastic,
  i.e., \\\sigma\_{u0}^2 = \exp(\gamma\_{u0}\[0\])\\.

- ln.var.v:

  optional one-sided formula; e.g. `ln.var.v = ~ z1 + z2`. Specifies
  exogenous variables entering the (log) variance of the random noise
  component. If `NULL`, the noise variance is homoskedastic, i.e.,
  \\\sigma\_{v0}^2 = \exp(\gamma\_{v0}\[0\])\\.

- skew.v:

  optional one-sided formula; e.g. `skew.v = ~ z5 + z6`. Allows the
  skewness of the noise to depend linearly on exogenous variables. If
  `NULL`, the skewness is constant across units.

- mean.u:

  optional one-sided formula; e.g. `mean.u = ~ z7 + z8`. Specifies
  whether the mean of the pre-truncated normal distribution of the
  inefficiency term is a linear function of exogenous variables. In
  cross-sectional models, used only when `distribution = "t"`. If
  `NULL`, the mean is constant across units. To be implemented.

- optim.report:

  integer. Verbosity level for reporting during optimization (if
  implemented). Default is `1`.

- optim.trace:

  integer. Trace level for optimization (if implemented). Default is
  `1`.

- reltol:

  numeric. Relative convergence tolerance used when maximizing the
  log-likelihood with `optim`. The algorithm stops if it cannot reduce
  the objective by a factor of `reltol * (abs(val) + reltol)` at a step.
  Default is `sqrt(.Machine$double.eps)`.

- lmtol:

  numeric. Convergence tolerance based on the scaled gradient (when
  applicable). Default is `1e-5`.

- maxit:

  numeric. Maximum number of iterations for the optimizer. Default is
  `199`.

- print.level:

  integer. Printing level: `1`—estimation results; `2`—optimization
  details; `3`—summary of (cost/technical) efficiencies;
  `4`—unit-specific point and interval estimates of efficiencies.
  Default is `3`.

- digits:

  integer. Number of digits for displaying estimates and efficiencies.
  Default is `4`.

- ...:

  additional arguments passed to internal methods or to `optim`, as
  relevant (e.g., `cost.eff.less.one = TRUE` for cost-frontier
  conventions).

- optim:

  logical. If `TRUE`, estimation proceeds via
  [`stats::optim`](https://rdrr.io/r/stats/optim.html); if `FALSE`, an
  internal routine (if provided) would be used. Default is `FALSE`.

- optim.method:

  character. Method passed to
  [`stats::optim`](https://rdrr.io/r/stats/optim.html) when
  `optim = TRUE`. Default is `"bfgs"`.

- optim.reltol:

  numeric. Relative tolerance specifically for `optim`; default `1e-8`.

## Value

An object of class `"snreg"` with maximum-likelihood estimates and
diagnostics:

- `par`:

  Numeric vector of ML parameter estimates at the optimum.

- `coef`:

  Named numeric vector equal to `par`.

- `vcov`:

  Variance–covariance matrix of the estimates.

- `sds`:

  Standard errors, `sqrt(diag(vcov))`.

- `ctab`:

  Coefficient table with columns `Coef.`, `SE`, `z`, `P>|z|`.

- `ll`:

  Maximized log-likelihood value.

- `hessian`:

  (When computed) Observed Hessian or OPG used to form `vcov`.

- `value`:

  (Optim-only, before aliasing) Objective value from `optim`.

- `counts`:

  (Optim-only) Iteration and evaluation counts from `optim`.

- `convergence`:

  Convergence code).

- `message`:

  (Optim-only) Message returned by `optim`, if any.

- `gradient`:

  (NR-only) Gradient at the solution.

- `gg`:

  (NR-only) Gradient-related diagnostic.

- `gHg`:

  (NR-only) Newton-step diagnostic.

- `theta_rel_ch`:

  (NR-only) Relative parameter change metric across iterations.

- `resid`:

  Regression residuals.

- `RSS`:

  Residual sum of squares `crossprod(resid)`.

- `shat2`:

  Residual variance estimate `var(resid)`.

- `shat`:

  Residual standard deviation `sqrt(shat2)`.

- `aic`:

  Akaike Information Criterion.

- `bic`:

  Bayesian Information Criterion.

- `Mallows`:

  Mallows' \\C_p\\-like statistic.

- `u`:

  Estimated inefficiency term (vector). Returned for models with an
  inefficiency component (e.g., exponential).

- `eff`:

  Efficiency scores `exp(-u)` (technical or cost, depending on `prod`).

- `sv`:

  Estimated (possibly unit-specific) standard deviation of the noise
  term.

- `su`:

  Estimated (possibly unit-specific) standard deviation or scale of the
  inefficiency term. For exponential models.

- `skewness`:

  Estimated skewness index (e.g., from the skewness equation).

- `esample`:

  Logical vector marking observations used in estimation.

- `n`:

  Number of observations used.

The returned object has class `"snreg"`.

## Details

Stochastic Frontier Model with a Skew-Normally Distributed Error Term

Models for `snsf` are specified symbolically. A typical model has the
form `y ~ x1 + ...`, where `y` represents the logarithm of outputs or
total costs and `{x1, ...}` is a set of inputs (for production) or
outputs and input prices (for cost), all typically in logs.

Options `ln.var.u` and `ln.var.v` allow for multiplicative
heteroskedasticity in the inefficiency and/or noise components; i.e.,
their variances can be modeled as exponential functions of exogenous
variables (including an intercept), as in Caudill et al. (1995).

## References

Badunenko, O., & Henderson, D. J. (2023). *Production analysis with
asymmetric noise*. Journal of Productivity Analysis, **61**(1), 1–18.
https://doi.org/10.1007/s11123-023-00680-5

Caudill, S. B., Ford, J. M., & Gropper, D. M. (1995). *Frontier
estimation and firm-specific inefficiency measures in the presence of
heteroskedasticity*. Journal of Business & Economic Statistics,
**13**(1), 105–111.

## See also

[`sf`](https://rdrr.io/pkg/npsf/man/sf.html)

## Author

Oleg Badunenko \<Oleg.Badunenko.brunel.ac.uk\>

## Examples

``` r
if (FALSE) { # \dontrun{

library(snreg)

data("banks07")
head(banks07)

myprod <- FALSE

# Translog cost function
spe.tl <- log(TC) ~ (log(Y1) + log(Y2) + log(W1) + log(W2))^2 +
  I(0.5 * log(Y1)^2) + I(0.5 * log(Y2)^2) +
  I(0.5 * log(W1)^2) + I(0.5 * log(W2)^2)


# -------------------------------------------------------------
# Specification 1: homoskedastic & symmetric
# -------------------------------------------------------------
formSV <- NULL   # variance equation
formSK <- NULL   # skewness equation
formSU <- NULL   # inefficiency equation (unused here)

m1 <- snsf(
  formula  = spe.tl,
  data     = banks07,
  prod     = myprod,
  ln.var.v = formSV,
  skew.v   = formSK
)

coef(m1)


# -------------------------------------------------------------
# Specification 2: heteroskedastic + skewed noise
# -------------------------------------------------------------
formSV <- ~ log(TA)      # heteroskedastic variance
formSK <- ~ ER           # skewness driver
formSU <- ~ LA + ER      # inefficiency

m2 <- snsf(
  formula  = spe.tl,
  data     = banks07,
  prod     = myprod,
  ln.var.v = formSV,
  skew.v   = formSK
)

coef(m2)

} # }
```
