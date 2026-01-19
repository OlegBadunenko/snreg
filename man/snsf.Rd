\name{snsf}
\alias{snsf}
\title{
 Stochastic Frontier Model with a Skew-Normally Distributed Error Term
}
\description{
 \code{snsf} performs maximum likelihood estimation of the parameters and technical or cost efficiencies Stochastic Frontier Model with a Skew-Normally Distributed Error Term.
}
\usage{
snsf(
  formula, data, subset, prod = TRUE,
  distribution,
  lmtol     = 1e-5,  maxit = 199, 
  reltol    = sqrt(.Machine$double.eps),
  start.val = NULL,
  init.sk   = 0.5,
  ln.var.u  = NULL,
  ln.var.v  = NULL,
  v.skew    = NULL,
  mean.u    = NULL,
  optim     = FALSE,
  optim.method = 'Nelder-Mead',
  optim.report = 1,
  optim.trace  = 1,
  optim.reltol = 1e-8,
  print.level = 3,
  ...)
}
\arguments{
 \item{formula}{
  an object of class ``formula'' (or one that can be coerced to that class): a symbolic description of the model. The details of model specification are given under `Details'.
 }
 \item{data}{
an optional data frame containing the variables in the model. If not found in data, the variables are taken from environment (\code{formula}), typically the environment from which \code{teradial} is called.
}
 \item{subset}{
an optional vector specifying a subset of observations for which technical efficiency is to be computed.
}
\item{prod}{
 logical. If \code{TRUE}, the estimates of parameters of stochastic production frontier model and of technical efficiencies are returned; if \code{FALSE}, the estimates of parameters of stochastic cost frontier model and of cost efficiencies are returned.
 }
 \item{distribution}{
 either \code{"h"} (half-normal), \code{"t"} (truncated normal), or \code{"e"} (exponential), specifying the distribution of inefficiency term.
 }
 \item{ln.var.v}{
 one-sided formula; e.g. \code{ln.var.v ~ z1 + z2}. Specifies exogenous variables entering the expression for variance of random effects component. If \code{NULL}, inefficiency term is assumed to be homoskedastic, i.e. \eqn{\sigma_v0^2 = exp(\gamma_v0[0])}.
 }
 \item{ln.var.u}{
 one-sided formula; e.g. \code{ln.var.u ~ z3 + z4}. Specifies exogenous variables entering the expression for variance of persistent inefficiency component. If \code{NULL}, inefficiency term is assumed to be homoskedastic, i.e. \eqn{\sigma_u0^2 = exp(\gamma_u0[0])}.
 }
 \item{v.skew}{
 one-sided formula; e.g. \code{v.skew ~ z5 + z6}. Specifies whether the skewness of the noise is a linear function of exogenous variables. If \code{NULL}, mean is assumed to be constant for all ids. 
 }
 \item{mean.u}{
 one-sided formula; e.g. \code{mean.u ~ z7 + z8}. Specifies whether the mean of pre-truncated normal distribution of inefficiency term is a linear function of exogenous variables. In cross-sectional models, used only when \code{distribution = "t"}. If \code{NULL}, mean is assumed to be constant for all ids. 
 }
\item{start.val}{
 numeric. Vector of starting values.
 }
  \item{lmtol}{
  numeric. Convergence tolerance based on scaled gradient. Defaults to \code{10^-5}.
 }
  \item{reltol}{
  numeric. Option for maximization of log likelihood using \code{optim}. Relative convergence tolerance. The algorithm stops if it is unable to reduce the value by a factor of \code{reltol * (abs(val) + reltol)} at a step. Defaults to \code{10^-15}.
 }
  \item{maxit}{
  numeric. The maximum number of iterations. Default is \code{3501}.
 }
 \item{digits}{
 numeric. Number of digits to be displayed in estimation results and for efficiency estimates. Default is 4.
 }
 \item{print.level}{
  numeric. 1 - print estimation results. 2 - print optimization details. 3 - print summary of point estimates of technical or cost efficiencies. 4 - print unit-specific point and interval estimates of technical or cost efficiencies. Default is 2.
 }
  }


\details{
 Models for \code{snsf} are specified symbolically. A typical model has the form \code{y ~ x1 + ...}, where \code{y} represents the logarithm of outputs or total costs and \code{\{x1,...\}} is a series of inputs or outputs and input prices (in logs).

Options \code{ln.var.u} and \code{ln.var.v} can be used if multiplicative heteroskedasticity of either inefficiency or random noise component (or both) is assumed; i.e. if their variances can be expressed as exponential functions of (e.g. size-related) exogenous variables (including intercept) (see Caudill et al. 1995).

}

\value{
 \code{snsf} returns a list of class \code{snsf} containing the following elements:
  \item{coef}{
   numeric. Named vector of ML parameter estimates.
  }
 \item{vcov}{
  matrix. Estimated covariance matrix of ML estimator.
 }
 \item{loglik}{
  numeric. Value of log-likelihood at ML estimates.
 }
 \item{efficiencies}{
  data frame. Contains point estimates of unit-specific technical or cost efficiencies:  exp(-E(u|e)) of Jondrow et al. (1982), E(exp(-u)|e) of Battese and Coelli (1988), and exp(-M(u|e)), where M(u|e) is the mode of conditional distribution of inefficiency term. In addition, estimated lower and upper bounds of (1-\eqn{\alpha})100\% two-sided prediction intervals are returned.
 }
 \item{su}{
  matrix. Estimated unit-specific variances of inefficiency term. Returned if \code{ln.var.u} is not \code{NULL}.
 }
  \item{sv}{
   matrix. Estimated unit-specific variances of random noise component. Returned if \code{ln.var.v} is not \code{NULL}.
 }
  \item{mu}{
   matrix. Estimated unit-specific means of pre-truncated normal distribution of inefficiency term. Returned if \code{tmean} is not \code{NULL}.
 }
 \item{esample}{
  logical. Returns TRUE if the observation in user supplied data is in the estimation subsample and FALSE otherwise.
 }
}


\references{

: @@@ (20@@), @title@. \emph{ @Journal@ }, \bold{@}(@), @−@. DOI: @@@

}

\author{
 Anonymouso@authors.com
}

\seealso{
\code{\link[npsf]{sf}}
}

\examples{
require( snsf )

data("banks07")
head(banks07)

\dontrun{

myprod <- FALSE

spe.tl <- log(TC) ~ (log(Y1) + log(Y2) + log(W1) + log(W2))^2 +
  I(0.5*log(Y1)^2) + I(0.5*log(Y2)^2) + I(0.5*log(W1)^2) + I(0.5*log(W2)^2) 

# Specification 1

formSV <- NULL
formSK <- NULL
formSU <- NULL

# Normal-Exponential model

m1.n <- npsf::sf  (
  formula = spe.tl, data = banks06_07, subset = year == 2007, prod = myprod,
  distribution = "e",
  ln.var.u.0i  = formSU,
  ln.var.v.0i  = formSV)

# Skewnormal-Exponential model

m1.sn <- snsf::snsf  (
  formula = spe.tl, data = banks06_07, subset = year == 2007, prod = myprod,
  ln.var.u  = formSU,
  ln.var.v  = formSV,
  v.skew    = formSK)

# Specification 2

formSK <- ~ ER

# Normal-Exponential model

m2.n <- npsf::sf  (
  formula = spe.tl, data = banks06_07, subset = year <= 2007, prod = myprod,
  distribution = "e",
  ln.var.u.0i  = formSU,
  ln.var.v.0i  = formSV,
  cost.eff.less.one = TRUE)

# Skewnormal-Exponential model

m2.sn <- snsf::snsf  (
  formula = spe.tl, data = banks06_07, subset = year <= 2007, prod = myprod,
  ln.var.u  = formSU,
  ln.var.v  = formSV,
  v.skew    = formSK,
  cost.eff.less.one = TRUE)

# Specification 3

formSV <- ~ log(TA)
formSK <- ~ ER
formSU <- ~ LA + ER

# Normal-Exponential model

m3.n <- npsf::sf  (
  formula = spe.tl, data = banks06_07, subset = year == 2007, prod = myprod,
  distribution = "e",
  ln.var.u.0i  = formSU,
  ln.var.v.0i  = formSV,
  cost.eff.less.one = TRUE)

# Skewnormal-Exponential model

m3.sn <- snsf::snsf  (
  formula = spe.tl, data = banks06_07, subset = year == 2007, prod = myprod,
  ln.var.u  = formSU,
  ln.var.v  = formSV,
  v.skew    = formSK)

}

}

\keyword{ Stochastic Frontier Analysis }
\keyword{ Heteroskedasticity }
\keyword{ Parametric efficiency analysis }

