# Linear Regression with Skew-Normal Errors

`snreg` fits a linear regression model where the disturbance term
follows a skew-normal distribution. The function supports multiplicative
heteroskedasticity of the noise variance via a log-linear specification
(`ln.var.v`) and allows the skewness parameter to vary linearly with
exogenous variables (`skew.v`).

## Usage

``` r
snreg(
  formula,
  data,
  subset,
  init.sk = NULL,
  ln.var.v = NULL,
  skew.v = NULL,
  start.val = NULL,
  technique = c("nr"),
  vcetype = c("aim"),
  lmtol = 1e-05,
  reltol = 1e-12,
  maxit = 199,
  optim.report = 1,
  optim.trace = 1,
  print.level = 3,
  digits = 4,
  only.data = FALSE,
  ...
)
```

## Arguments

- formula:

  an object of class `formula` specifying the regression: typically
  `y ~ x1 + ...`, where `y` is the dependent variable and `x`'s are
  regressors.

- data:

  an optional `data.frame` containing the variables in `formula`. If not
  found in `data`, variables are taken from `environment(formula)`.

- subset:

  an optional logical or numeric vector specifying the subset of
  observations to be used in estimation.

- init.sk:

  numeric. Initial value for the (global) skewness parameter of the
  noise; can be `NULL` if `skew.v` is supplied with its own coefficients
  to initialize.

- ln.var.v:

  optional one-sided formula; e.g. `ln.var.v ~ z1 + z2`. Specifies
  exogenous variables entering the (log) variance of the random noise
  component. If `NULL`, the noise variance is homoskedastic.

- skew.v:

  optional one-sided formula; e.g. `skew.v ~ z3 + z4`. Specifies
  exogenous variables determining the skewness of the noise via a linear
  index; if `NULL`, the skewness is constant (scalar).

- start.val:

  optional numeric vector of starting values for all free parameters
  (regression coefficients, variance/heteroskedasticity parameters,
  skewness parameters).

- technique:

  character vector giving the preferred maximization routine(s) in order
  of preference. Currently recognized keywords include `"nr"`
  (Newton–Raphson), `"bhhh"`, `"nm"` (Nelder–Mead), `"bfgs"`, `"cg"`.
  This scaffold does not implement them yet, but records the choice.

- vcetype:

  character specifying the variance-covariance estimator type: `"aim"`
  for the approximated information matrix or `"opg"` for the outer
  product of gradients. Default is `"aim"`.

- lmtol:

  numeric. Convergence tolerance based on the scaled gradient (if
  applicable). Default is `1e-5`.

- reltol:

  numeric. Relative convergence tolerance for likelihood maximization.
  Default is `1e-12`.

- maxit:

  integer. Maximum number of iterations for the optimizer. Default is
  `199`.

- optim.report:

  integer. Verbosity for reporting progress (if implemented). Default is
  `1`.

- optim.trace:

  integer. If positive, tracing information is printed (if implemented).
  Default is `1`.

- print.level:

  integer. Printing level for summaries: `1`—print estimation results;
  `2`—print optimization details; `3`—print compact summary. Default
  `3`.

- digits:

  integer. Number of digits for printing. Default `4`.

- only.data:

  logical. If `TRUE`, the function returns only the constructed model
  matrices and design sets (no estimation). Default `FALSE`.

- ...:

  additional arguments reserved for future methods (e.g., box
  constraints).

## Value

An object of class `"snreg"` containing the maximum-likelihood results
and, depending on the optimization routine, additional diagnostics:

- `par`:

  Numeric vector of parameter estimates at the optimum.

- `coef`:

  Named numeric vector equal to `par`.

- `vcov`:

  Variance–covariance matrix of the estimates.

- `sds`:

  Standard errors, computed as `sqrt(diag(vcov))`.

- `ctab`:

  Coefficient table with columns: `Estimate`, `Std.Err`, `Z value`,
  `Pr(>z)`.

- `RSS`:

  Residual sum of squares.

- `esample`:

  Logical vector indicating which observations were used in estimation.

- `n`:

  Number of observations used in the estimation sample.

- `skewness`:

  Vector of the fitted skewness index.

- `hessian`:

  (BFGS only) Observed Hessian at the optimum. If `vcetype == "opg"`,
  this is set to the negative outer product of the individual gradients;
  otherwise a numerical Hessian is computed.

- `value`:

  (BFGS only) Objective value returned by `optim`. With
  `control$fnscale = -1`, this equals the maximized log-likelihood.

- `counts`:

  (BFGS only) Number of iterations / function evaluations returned by
  `optim`.

- `convergence`:

  (BFGS only) Convergence code from `optim`.

- `message`:

  (BFGS only) Additional `optim` message, if any.

- `ll`:

  Maximized log-likelihood value.

- `gradient`:

  (NR only) Gradient at the solution.

- `gg`:

  (NR only) Optional gradient-related diagnostic.

- `gHg`:

  (NR only) Optional Newton-step diagnostic.

- `theta_rel_ch`:

  (NR only) Relative parameter change metric across iterations.

The returned object has class `"snreg"`.

## Details

Linear Regression with Skew-Normal Errors

The model is \$\$y_i = x_i^\top \beta + \varepsilon_i,\quad
\varepsilon_i \sim SN(0, \sigma_i^2, \alpha_i),\$\$ where \\SN\\ denotes
the skew-normal distribution (Azzalini).

Heteroskedasticity in the noise variance (if specified via `ln.var.v`)
is modeled as \$\$\log(\sigma_i^2) = w_i^\top \gamma_v,\$\$ and the
(optional) covariate-driven skewness (if specified via `skew.v`) as
\$\$\alpha_i = s_i^\top \delta.\$\$

This function constructs the model frame and design matrices for
\\\beta\\, \\\gamma_v\\, and \\\delta\\, and is designed to be paired
with a maximum likelihood routine to estimate parameters and
(optionally) their asymptotic covariance via either AIM or OPG.

## References

Azzalini, A. (1985). *A Class of Distributions Which Includes the Normal
Ones*. Scandinavian Journal of Statistics, 12(2), 171–178.

Azzalini, A., & Capitanio, A. (2014). *The Skew-Normal and Related
Families*. Cambridge University Press.

## Examples

``` r
if (FALSE) { # \dontrun{

library(snreg)

data("banks07")
head(banks07)

# Translog cost function
spe.tl <- log(TC) ~ (log(Y1) + log(Y2) + log(W1) + log(W2))^2 +
  I(0.5 * log(Y1)^2) + I(0.5 * log(Y2)^2) +
  I(0.5 * log(W1)^2) + I(0.5 * log(W2)^2)


# -------------------------------------------------------------
# Specification 1: homoskedastic & symmetric noise
# -------------------------------------------------------------
formSV <- NULL     # variance equation
formSK <- NULL     # skewness equation

m1 <- snreg(
  formula  = spe.tl,
  data     = banks07,
  ln.var.v = formSV,
  skew.v   = formSK
)

coef(m1)


# -------------------------------------------------------------
# Specification 2: heteroskedastic + skewed noise
# -------------------------------------------------------------
formSV <- ~ log(TA)   # heteroskedasticity in v
formSK <- ~ ER        # skewness driven by equity ratio

m2 <- snreg(
  formula  = spe.tl,
  data     = banks07,
  ln.var.v = formSV,
  skew.v   = formSK
)

coef(m2)

} # }
```
