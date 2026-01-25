
#' Stochastic Frontier Model with a Skew-Normally Distributed Error Term
#'
#' @title Stochastic Frontier Model with a Skew-Normally Distributed Error Term
#'
#' @description
#' \code{snsf} performs maximum likelihood estimation of the parameters and
#' technical or cost efficiencies in a Stochastic Frontier Model with a
#' skew-normally distributed error term.
#'
#' @param formula
#' an object of class \code{formula} specifying the frontier:
#' a typical model is \code{y ~ x1 + ...}, where \code{y} is the
#' log of output (or total cost), and \code{x}'s are inputs (or outputs and
#' input prices, in logs). See \strong{Details}.
#'
#' @param data
#' an optional \code{data.frame} containing the variables in \code{formula}.
#' If not found in \code{data}, variables are taken from \code{environment(formula)}.
#'
#' @param subset
#' an optional logical or numeric vector specifying a subset of observations
#' for which the model is estimated and efficiencies are computed.
#'
#' @param prod
#' logical. If \code{TRUE}, estimates correspond to a stochastic \emph{production}
#' frontier and technical efficiencies are returned; if \code{FALSE}, estimates
#' correspond to a stochastic \emph{cost} frontier and cost efficiencies are returned.
#' Default is \code{TRUE}.
#'
#' @param distribution
#' character scalar specifying the distribution of the inefficiency term: default \code{"e"} (exponential).
#' \code{"h"} (half-normal) and \code{"t"} (truncated normal) to be implemented.
#'
#' @param lmtol
#' numeric. Convergence tolerance based on the scaled gradient (when applicable).
#' Default is \code{1e-5}.
#'
#' @param maxit
#' numeric. Maximum number of iterations for the optimizer. Default is \code{199}.
#'
#' @param reltol
#' numeric. Relative convergence tolerance used when maximizing the log-likelihood
#' with \code{optim}. The algorithm stops if it cannot reduce the objective by
#' a factor of \code{reltol * (abs(val) + reltol)} at a step. Default is
#' \code{sqrt(.Machine$double.eps)}.
#'
#' @param start.val
#' optional numeric vector of starting values for the optimizer.
#'
#' @param init.sk
#' numeric. Initial value for the skewness parameter of the noise component;
#' default is \code{0.5}.
#'
#' @param ln.var.u
#' optional one-sided formula; e.g. \code{ln.var.u = ~ z3 + z4}. Specifies exogenous
#' variables entering the (log) variance of the inefficiency component. If
#' \code{NULL}, the inefficiency variance is homoskedastic, i.e.,
#' \eqn{\sigma_{u0}^2 = \exp(\gamma_{u0}[0])}.
#'
#' @param ln.var.v
#' optional one-sided formula; e.g. \code{ln.var.v = ~ z1 + z2}. Specifies exogenous
#' variables entering the (log) variance of the random noise component. If
#' \code{NULL}, the noise variance is homoskedastic, i.e.,
#' \eqn{\sigma_{v0}^2 = \exp(\gamma_{v0}[0])}.
#'
#' @param skew.v
#' optional one-sided formula; e.g. \code{skew.v = ~ z5 + z6}. Allows the skewness
#' of the noise to depend linearly on exogenous variables. If \code{NULL}, the
#' skewness is constant across units.
#'
#' @param mean.u
#' optional one-sided formula; e.g. \code{mean.u = ~ z7 + z8}. Specifies whether the
#' mean of the pre-truncated normal distribution of the inefficiency term is a
#' linear function of exogenous variables. In cross-sectional models, used only
#' when \code{distribution = "t"}. If \code{NULL}, the mean is constant across units. To be implemented.
#'
#' @param optim
#' logical. If \code{TRUE}, estimation proceeds via \code{stats::optim}; if
#' \code{FALSE}, an internal routine (if provided) would be used. Default is \code{FALSE}.
#'
#' @param optim.method
#' character. Method passed to \code{stats::optim} when \code{optim = TRUE}.
#' Default is \code{"bfgs"}.
#'
#' @param optim.report
#' integer. Verbosity level for reporting during optimization (if implemented).
#' Default is \code{1}.
#'
#' @param optim.trace
#' integer. Trace level for optimization (if implemented). Default is \code{1}.
#'
#' @param optim.reltol
#' numeric. Relative tolerance specifically for \code{optim}; default \code{1e-8}.
#'
#' @param print.level
#' integer. Printing level: \code{1}—estimation results;
#' \code{2}—optimization details; \code{3}—summary of (cost/technical)
#' efficiencies; \code{4}—unit-specific point and interval estimates of
#' efficiencies. Default is \code{3}.
#'
#' @param digits
#' integer. Number of digits for displaying estimates and efficiencies. Default is \code{4}.
#'
#' @param ...
#' additional arguments passed to internal methods or to \code{optim}, as relevant
#' (e.g., \code{cost.eff.less.one = TRUE} for cost-frontier conventions).
#'
#' @details
#' Models for \code{snsf} are specified symbolically. A typical model has the form
#' \code{y ~ x1 + ...}, where \code{y} represents the logarithm of outputs or total
#' costs and \code{\{x1, ...\}} is a set of inputs (for production) or outputs and
#' input prices (for cost), all typically in logs.
#'
#' Options \code{ln.var.u} and \code{ln.var.v} allow for multiplicative
#' heteroskedasticity in the inefficiency and/or noise components; i.e., their
#' variances can be modeled as exponential functions of exogenous variables
#' (including an intercept), as in Caudill et al. (1995).
#'

#' @return
#' An object of class \code{"snreg"} with maximum-likelihood estimates and diagnostics:
#'
#' \describe{
#'
#'   \item{\code{par}}{Numeric vector of ML parameter estimates at the optimum.}
#'   \item{\code{coef}}{Named numeric vector equal to \code{par}.}
#'   \item{\code{vcov}}{Variance–covariance matrix of the estimates.}
#'   \item{\code{sds}}{Standard errors, \code{sqrt(diag(vcov))}.}
#'   \item{\code{ctab}}{Coefficient table with columns
#'         \code{Coef.}, \code{SE}, \code{z}, \code{P>|z|}.}
#'
#'   \item{\code{ll}}{Maximized log-likelihood value.}
#'   \item{\code{hessian}}{(When computed) Observed Hessian or OPG used to form \code{vcov}.}
#'   \item{\code{value}}{(Optim-only, before aliasing) Objective value from \code{optim}.}
#'   \item{\code{counts}}{(Optim-only) Iteration and evaluation counts from \code{optim}.}
#'   \item{\code{convergence}}{Convergence code).}
#'   \item{\code{message}}{(Optim-only) Message returned by \code{optim}, if any.}
#'   \item{\code{gradient}}{(NR-only) Gradient at the solution.}
#'   \item{\code{gg}}{(NR-only) Gradient-related diagnostic.}
#'   \item{\code{gHg}}{(NR-only) Newton-step diagnostic.}
#'   \item{\code{theta_rel_ch}}{(NR-only) Relative parameter change metric across iterations.}
#'
#'   \item{\code{resid}}{Regression residuals.}
#'   \item{\code{RSS}}{Residual sum of squares \code{crossprod(resid)}.}
#'   \item{\code{shat2}}{Residual variance estimate \code{var(resid)}.}
#'   \item{\code{shat}}{Residual standard deviation \code{sqrt(shat2)}.}
#'   \item{\code{aic}}{Akaike Information Criterion.}
#'   \item{\code{bic}}{Bayesian Information Criterion.}
#'   \item{\code{Mallows}}{Mallows' \eqn{C_p}-like statistic.}
#'
#'   \item{\code{u}}{Estimated inefficiency term (vector). Returned for models with
#'         an inefficiency component (e.g., exponential).}
#'   \item{\code{eff}}{Efficiency scores \code{exp(-u)} (technical or cost, depending on \code{prod}).}
#'   \item{\code{sv}}{Estimated (possibly unit-specific) standard deviation of the noise term.}
#'   \item{\code{su}}{Estimated (possibly unit-specific) standard deviation or scale of the
#'         inefficiency term. For exponential models.}
#'   \item{\code{skewness}}{Estimated skewness index (e.g., from the skewness equation).}
#'
#'   \item{\code{esample}}{Logical vector marking observations used in estimation.}
#'   \item{\code{n}}{Number of observations used.}
#' }
#'
#' The returned object has class \code{"snreg"}.
#'
#' @references
#' Badunenko, O., & Henderson, D. J. (2023).
#' \emph{Production analysis with asymmetric noise}.
#' Journal of Productivity Analysis, \bold{61}(1), 1–18.
#' https://doi.org/10.1007/s11123-023-00680-5
#'
#' Caudill, S. B., Ford, J. M., & Gropper, D. M. (1995).
#' \emph{Frontier estimation and firm-specific inefficiency measures in the presence of heteroskedasticity}.
#' Journal of Business & Economic Statistics, \bold{13}(1), 105–111.
#'
#' @author
#' Oleg Badunenko <Oleg.Badunenko.brunel.ac.uk>
#'
#' @seealso \code{\link[npsf]{sf}}
#'

#' @examples
#' \dontrun{
#'
#' library(snreg)
#'
#' data("banks07")
#' head(banks07)
#'
#' myprod <- FALSE
#'
#' # Translog cost function
#' spe.tl <- log(TC) ~ (log(Y1) + log(Y2) + log(W1) + log(W2))^2 +
#'   I(0.5 * log(Y1)^2) + I(0.5 * log(Y2)^2) +
#'   I(0.5 * log(W1)^2) + I(0.5 * log(W2)^2)
#'
#'
#' # -------------------------------------------------------------
#' # Specification 1: homoskedastic & symmetric
#' # -------------------------------------------------------------
#' formSV <- NULL   # variance equation
#' formSK <- NULL   # skewness equation
#' formSU <- NULL   # inefficiency equation (unused here)
#'
#' m1 <- snsf(
#'   formula  = spe.tl,
#'   data     = banks07,
#'   prod     = myprod,
#'   ln.var.v = formSV,
#'   skew.v   = formSK
#' )
#'
#' coef(m1)
#'
#'
#' # -------------------------------------------------------------
#' # Specification 2: heteroskedastic + skewed noise
#' # -------------------------------------------------------------
#' formSV <- ~ log(TA)      # heteroskedastic variance
#' formSK <- ~ ER           # skewness driver
#' formSU <- ~ LA + ER      # inefficiency (not used by snsf if prod=FALSE)
#'
#' m2 <- snsf(
#'   formula  = spe.tl,
#'   data     = banks07,
#'   prod     = myprod,
#'   ln.var.v = formSV,
#'   skew.v   = formSK
#' )
#'
#' coef(m2)
#'
#' }
#'
#' @keywords Stochastic Frontier Analysis Heteroskedasticity Parametric efficiency analysis
#' @export
snsf <- function(
    formula, data, subset,
    distribution = "e",
    prod      = TRUE,
    start.val = NULL,
    init.sk   = NULL,#.5*(2*prod-1),
    ln.var.u  = NULL,
    ln.var.v  = NULL,
    skew.v    = NULL,
    mean.u    = NULL,
    technique = c('nr'),#,'bhhh','nm', 'bfgs', 'cg'),
    vcetype   = c('aim'),#, 'opg'), # `approximated information matrix` or `outer product of gradients`
    optim.report = 1,
    optim.trace  = 1,
    reltol = 1e-12,
    lmtol = 1e-5,  
    maxit = 199,
    print.level = 3, 
    report    = 1,
    trace     = 1,
    threads   = 1,
    only.data = FALSE,
    digits = 4,
    ...
) {

  # threads <- 1
  
  # cat("0\n")
  myprod <- ifelse(prod, 1, -1)
  
  # handle distribution
  
  if(length(distribution) != 1){
    stop("Distribution of inefficiency term should be specified.")
  } else {
    distribution <- tolower(substr(distribution, 1, 1 ))
  }
  
  if( !distribution %in% c("t","h","e") ){
    stop("'distribution' is invalid")
  }
  
  # handle technique
  
  technique <- technique[1]
  
  if( !technique %in% c('nr','bhhh','nm', 'bfgs','cg') ){
    stop("'technique' is invalid")
  }
  
  if( !vcetype %in% c('aim', 'opg') ){
    stop("'vcetype' is invalid")
  }
  
  if(technique == 'nm') technique <- 'Nelder-Mead'
  if(technique == 'bfgs') technique <- 'BFGS'
  if(technique == 'cg') technique <- 'CG'
  
  
  # if(distribution == "h"){
  #   mean.u.zero  = TRUE
  # } else {
  #   mean.u.zero  = FALSE
  # }
  
  if(is.null(mean.u) == FALSE & distribution != "t"){
    stop("Option 'mean.u' can be used only when distribution of inefficiency term is truncated normal.")
  }
  
  # * data --------------------------------------------------------------------
  
  # # cat.print(data[subset,])
  # yxz <- sn.tn.get.data(formX=formula, data=data, subset=subset, formNV=ln.var.v, formNS=skew.v, formIV=ln.var.u, formIL=mean.u)
  # cat("1\n")
  
  mf0 <- match.call(expand.dots = FALSE, call = sys.call(sys.parent(n = 0)))
  if(is.null(ln.var.v)){
    ln.var.v <- ~ 1
    ksv <- 1
    # cat.print(ksv)
  } else {
    ksv <- 17
  }
  if(is.null(skew.v)){
    skew.v <- ~ 1
    ksk <- 1
    # cat.print(ksk)
  } else {
    ksk <- 17
  }
  if(is.null(ln.var.u)){
    ln.var.u <- ~ 1
    ksu <- 1
    # cat.print(ksu)
  } else {
    ksu <- 17
  }
  if(distribution == "t"){
    if(is.null(mean.u)){
      mean.u <- ~ 1
      kmu <- 1
      # cat.print(kmu)
    }else {
      kmu <- 17
    }
  }
  # cat("03\n")
  if(distribution == "t"){
    form1 <- as.Formula(formula, ln.var.v, skew.v, ln.var.u, mean.u)
  } else {
    form1 <- as.Formula(formula, ln.var.v, skew.v, ln.var.u)
  }
  # form1 <- as.Formula(formula, ln.var.v, skew.v, ln.var.u)
  # cat.print(form1)
  # cat("04\n")
  
  data.order <- as.numeric(rownames(data)) # seq_len(nrow(data))
  
  mf <- mf0
  mf$formula <- form1 #formula( form )
  # cat("05\n")
  m <- match(c("formula", "data", "subset"), names(mf), 0L)
  # cat("06\n")
  # cat.print(m)
  mf <- mf[c(1L, m)]
  # cat("07\n")
  # cat.print(mf)
  mf[[1L]] <- as.name("model.frame")
  # cat("08\n")
  # cat.print(mf)
  # cat.print(needed.frame)
  # mf <- eval(mf, sys.frame(sys.parent(n = needed.frame+2)))
  mf <- eval(mf, parent.frame())
  # cat.print(mf)
  # cat("09\n")
  # cat.print(mf)
  mt <- attr(mf, "terms")
  x <- as.matrix(model.matrix(mt, mf))
  # print(x)
  # now get the names in the entire data
  # esample <- seq_len( nrow(data) ) %in% as.numeric(rownames(x))
  esample <- rownames(data) %in% rownames(x)
  # cat("10\n")
  # cat.print(mt)
  X <- model.matrix(mt, mf)
  # cat("11\n")
  # cat.print(X)
  esample.nu <- as.numeric(rownames(X))
  esample <- data.order %in% esample.nu
  
  Y <- as.matrix(model.part(form1, data = mf, lhs = 1, drop = FALSE))
  # print(length(Y))
  # print(Y)
  X <- as.matrix( model.matrix(formula(form1, lhs = 0, rhs = 1), data = mf))
  # print(X)
  n <- nrow(Y)
  k <- ncol(X)
  # ksv:  noise, scale/variance
  # zsv
  if(ksv == 1){
    zsv <- matrix(1, nrow = sum(esample), ncol = 1)
  }
  else {
    zsv <- as.matrix( model.matrix(formula(form1, lhs = 0, rhs = 2), data = mf))
    ksv <- ncol(zsv)
  }
  # zsk
  if(ksk == 1){
    zsk <- matrix(1, nrow = sum(esample), ncol = 1)
  }
  else {
    zsk <- as.matrix( model.matrix(formula(form1, lhs = 0, rhs = 3), data = mf))
    ksk <- ncol(zsk)
  }
  # zsu
  if(ksu == 1){
    zsu <- matrix(1, nrow = sum(esample), ncol = 1)
  }
  else {
    zsu <- as.matrix( model.matrix(formula(form1, lhs = 0, rhs = 4), data = mf))
    ksu <- ncol(zsu)
  }
  if(distribution == "t"){
    # zmu
    if(kmu == 1){
      zmu <- matrix(1, nrow = sum(esample), ncol = 1)
    }
    else {
      zmu <- as.matrix( model.matrix(formula(form1, lhs = 0, rhs = 5), data = mf))
      kmu <- ncol(zmu)
    }
  } else {
    kmu <- 0
  }
  
  # from `sf`
  
  colnames(X)[1]   <- "(Intercept)"
  colnames(zsv)[1] <- "(Intercept)"
  colnames(zsk)[1] <- "(Intercept)"
  colnames(zsu)[1] <- "(Intercept)"
  if(distribution == "t"){
    colnames(zmu)[1] <- "(Intercept)"
  }
  
  abbr.length <- 34
  names_x <- abbreviate(colnames(X), abbr.length+6, strict = TRUE, dot = FALSE, method = "both.sides")
  # names_zu <-  abbreviate(colnames(Zu), 9, strict = TRUE, dot = FALSE)
  # names_zv <-  abbreviate(colnames(Zv), 9, strict = TRUE, dot = FALSE)
  # Zu_colnames <- abbreviate(colnames(Zu), abbr.length, strict = TRUE, dot = FALSE, method = "both.sides")
  names_zu <- paste0("lnVARu0i_", abbreviate(colnames(zsu), abbr.length, strict = TRUE, dot = FALSE, method = "both.sides"))
  names_zv <- paste0("lnVARv0i_", abbreviate(colnames(zsv), abbr.length, strict = TRUE, dot = FALSE, method = "both.sides"))
  names_zk <- paste0("Skew_v0i_", abbreviate(colnames(zsk), abbr.length, strict = TRUE, dot = FALSE, method = "both.sides"))
  
  # cat.print(names_zu)
  # cat.print(names_zv)  
  # cat.print(names_zk)
  
  if(distribution == "t"){
    names_zm <- paste0("mean_u0i_", abbreviate(colnames(zmu), 9, strict = TRUE, dot = FALSE))
  }
  
  coef.names.full <- c(
    names_x,
    paste("lnVARv0i_",c("(Intercept)", names_zv[-1]),"", sep = ""),
    paste("Skew_v0i_",c("(Intercept)", names_zk[-1]),"", sep = ""),
    paste("lnVARu0i_",c("(Intercept)", names_zu[-1]),"", sep = "")
  )
  if(distribution == "t"){
    coef.names.full <- c(coef.names.full, paste("mean_u0i_",c("(Intercept)", names_zm[-1]),"", sep = ""))
  }
  # cat.print(coef.names.full)
  coef.names.full <- c(
    names_x,  names_zv, names_zk, names_zu
  )
  if(distribution == "t"){
    coef.names.full <- c(coef.names.full, names_zm)
  }
  # cat.print(coef.names.full)
  
  # return items
  # colnames(X)[1]   <- "(Intercept)"
  # colnames(zsv)[1] <- "(Intercept)_ln_VAR_vi"
  # colnames(zsk)[1] <- "(Intercept)_SKEW_vi"
  # colnames(zsu)[1] <- "(Intercept)_ln_VAR_ui"
  # if(distribution == "t"){
  #   colnames(zmu)[1] <- "(Intercept)_Umu"
  # }
  
  # * starting values ---------------------------------------------------------
  
  # Run a N-SN regression 0
  # if(ksv==1){
  #   theta0 <- c( coef(lm(Y~X-1)), log(.1) )
  # } else {
  #   theta0 <- c( coef(lm(Y~X-1)), log(.1), rep(0, ksv-1) )
  # }
  # 
  # if(ksk==1){
  #   theta0 <- c( theta0, -.1)
  # } else {
  #   theta0 <- c( theta0, -.1, rep(0, ksk-1))
  # }
  # names(theta0) <- coef.names.full[1:(k+ksv+ksk)]
  # 
  # # cat.print(theta0)
  # tymch1 <- optim(
  #   par = theta0, fn = .ll.sn., y = Y, x = X, zsv = zsv, zsk = zsk, k = k, ksv = ksv, ksk = ksk,
  #   method = "BFGS", control = list(fnscale = -1, trace = 1, maxit = 10000), 
  #   hessian = TRUE)
  
  tymch2 <- lm(Y ~ X - 1)
  olsres  <- resid(tymch2)
  shat <- mean(olsres^2)
  
  gamma1.0 <- min(.99, mean(olsres^3)/shat )
  
  r  <- sign(gamma1.0)*(2*abs(gamma1.0)/(4-pi))^(1/3)
  delta <- r/(sqrt(2/pi) *sqrt(1+r^2))
  if(abs(delta) > .99) delta <- sign(delta) * .99
  # print(delta)
  alpha.0 <- delta/sqrt(1-delta^2)
  
  theta0 <- c( coef(lm(Y~X-1)), `lnVARv0i_(Intercept)` = log(.1), `Skew_v0i_(Intercept)` = alpha.0)
  
  if(print.level >= 1){
    max.name.length <- max(nchar(names(theta0) ))
    est.rez.left <- floor( (max.name.length+42-33) / 2 )
    est.rez.right <- max.name.length+42-33 - est.rez.left
    cat("\n",rep("-", est.rez.left)," Regression with skewed errors: ",rep("-", est.rez.right),"\n\n", sep ="")
  }
  # mytrace <- ifelse(print.level < 2, 0, 1)
  if(print.level <= 2){
    trace <- trace1 <- 0
  } else {
    trace1 <- optim.trace
  }
  tymch1 <- optim(
    par = theta0, fn = .ll.sn0, y = Y, x = X, zsv = zsv, zsk = zsk, k = k, ksv = ksv, ksk = ksk,
    method = "BFGS", control = list(fnscale = -1, trace = trace1, REPORT = optim.report, maxit = 10000), 
    hessian = TRUE)
  
  # cat.print(tymch1)
  # table with coefs
  tymch1$par  <- c(tymch1$par, sv = unname(sqrt(exp(tymch1$par[k+1]))) )
  tymch1$sds  <- suppressWarnings( sqrt(diag(solve(-tymch1$hessian))) )
  tymch1$sds  <- c(tymch1$sds, tymch1$sds[k+1]*exp(tymch1$par[k+1]))
  
  tymch1$ctab <- cbind(tymch1$par, tymch1$sds, tymch1$par/tymch1$sds,2*(1-pnorm(abs(-tymch1$par/tymch1$sds))))
  colnames(tymch1$ctab) <- c("Estimate", "Std.Err", "Z value", "Pr(>z)")
  
  if(print.level >= 1){
    printCoefmat(tymch1$ctab)
    cat("",rep("_", max.name.length+42-1),"", "\n", sep = "")
  }
  # cat.print(tymch1)
  # return(tymch1)
  
  # Run a N-SN regression 1
  
  # Initial value of the skewness parameter
  
  if(abs(tymch1$par[k+2]) < .01) tymch1$par[k+2] <- .01 * sign(tymch1$par[k+2])
  
  tymch1$par[k+2] <- alpha.0
  
  if(is.null(init.sk)){
    # sk.initial <- tymch1$par[(k+ksv+1):(k+ksv+ksk)] * .5
    if(ksk==1){
      sk.initial <- tymch1$par[k+2] * .5
    } else {
      sk.initial <- c(tymch1$par[k+2] * .5, rep(0,ksk-1))
    }
  } else {
    if(ksk==1){
      sk.initial <- init.sk
    } else {
      sk.initial <- c(init.sk, rep(0,ksk-1))
    }
  }
  
  # cat.print(tymch1$par[(1):(k)])
  
  # /* Use MOM to obtain starting values */
  
  tymch2 <- lm(Y ~ X - 1)
  # summary(tymch1)
  theta0 <- coef(tymch2)
  theta0 <- tymch1$par[(1):(k)]
  # cat.print(theta0)
  # hist(resid(tymch1))
  olsres  <- resid(tymch2)
  olsres2 <- olsres^2; m2 = mean(olsres2)
  # cat.print(m2)
  olsres3 <- olsres^3; m3 = mean(olsres3)
  # cat.print(m3)
  # m3t     <- m3/sqrt(6*m2^3/n)
  # /* negative skewness, so one-side test */
  # p_m3t   <- pnorm(myprod*m3t)
  # /* correct 3rd moment to be negative */
  # cat.print(m3)
  m3      <- ifelse(m3*myprod < 0, m3, -0.0001*myprod)
  # cat.print(m3)
  
  if(distribution == "e"){
    ou2 <- (-myprod*m3/2)^(2/3)
    # cat.print(ou2)
    lnou2 <- log(ou2)
    # cat.print(lnou2)
    ov2 <- m2 - ou2
    # cat.print(ov2)
    lnov2 <- ifelse(ov2>0, log(ov2), log(.0001))
    # cat.print(lnov2)
    # theta0[1] <- theta0[1] + myprod*(sqrt(2/pi))
    # cat.print(theta0)
  } else {
    ou2 <- (myprod*m3/(sqrt(2/pi) *(1-4/pi)))^(2/3)
    lnou2 <- log(ou2)
    ov2 <- m2 - (1-2/pi) * ou2
    lnov2 <- ifelse(ov2>0, log(ov2), log(.0001))
    # theta0[1] <- theta0[1] + myprod*(sqrt(2/pi))*sqrt(ou2)
  }
  
  # some most probably bad sv
  # log SV2
  if(ksv == 1){
    theta0 <- c(theta0, lnov2)
  } else {
    theta0 <- c(theta0, coef(lm( rep(lnov2,n) ~ zsv - 1)))
  }
  # some most probably bad sk
  # if(ksk == 1){
  #   # theta0 <- c(theta0, myprod*init.sk)
  #   theta0 <- c(theta0, sk.initial)
  # } else {
  #   # theta0 <- c(theta0, myprod*init.sk, rep(0.0, (ksk-1) ))
  #   theta0 <- c(theta0, sk.initial, rep(0.0, (ksk-1) ))
  # }
  theta0 <- c(theta0, sk.initial)
  # some most probably bad su
  # log SU2
  if(ksu == 1){
    theta0 <- c(theta0, lnou2)
  } else {
    theta0 <- c(theta0, coef(lm( rep(lnou2,n) ~ zsu - 1)))
  }
  if(distribution == "t"){
    # some most probably bad mu
    if(kmu == 1){
      theta0 <- c(theta0, 0.01)
    } else {
      theta0 <- c(theta0, 0.01, rep(0, kmu-1 ))
    }
  }
  
  
  if(distribution == "t"){
    theta_names <- c(colnames(X), colnames(zsv), colnames(zsk), colnames(zsu), colnames(zmu))
    tn <- 1
  } else {
    theta_names <- c(colnames(X), colnames(zsv), colnames(zsk), colnames(zsu))
    tn <- 0
  }
  
  k.all <- length(theta_names)
  
  # cat.print(muIsZero)
  
  # cat.print(theta0)
  # cat.print(coef.names.full)
  names(theta0) <- coef.names.full
  
  if(!is.null(start.val)){
    names.start.val.2.use <- intersect( names(theta0), names(start.val) )
    # theta0[names(theta0) %in% names(start.val)] <- start.val
    theta0[names(theta0) %in% names.start.val.2.use] <- start.val[names.start.val.2.use]
  }
  
  if(print.level > 2){
    cat.print(theta0)
  }
  
  if(only.data){
    if(distribution == "t"){
      return(
        list(
          theta0 = theta0, myprod = myprod,
          y=Y, x=X, zsv=zsv, zsk=zsk, zsu=zsu, zmu=zmu,
          n.ids=n, k=k, ksv=ksv, ksk=ksk, ksu=ksu, kmu=kmu
        )
      )
    } else {
      return(
        list(
          theta0 = theta0, myprod = myprod,
          y=Y, x=X, zsv=zsv, zsk=zsk, zsu=zsu,
          n.ids=n, k=k, ksv=ksv, ksk=ksk, ksu=ksu
        )
      )
    }
  }
  
  # ll0 <- .ll.sn.exp(theta0, myprod = myprod, 
  #                   y=yxz$y, x=yxz$x, zsv=yxz$zsv, zsk=yxz$zsk, zsu=yxz$zsu, 
  #                   n.ids=yxz$n, k=yxz$k, ksv=yxz$ksv, ksk=yxz$ksk, ksu=yxz$ksu )
  
  # cat.print(ll0)
  
  # g.appro <- .gr.sn.exp.appr(theta0, myprod=myprod, y=Y, x=X, zsv=zsv, zsk=zsk, zsu=zsu, k=k, ksv=ksv, ksk=ksk, ksu=ksu, n.ids=n)
  # g.analy <- .gr.sn.exp(theta0, myprod=myprod, y=Y, x=X, zsv=zsv, zsk=zsk, zsu=zsu, k=k, ksv=ksv, ksk=ksk, ksu=ksu, n.ids=n)
  # 
  # cat.print(g.appro)
  # cat.print(g.analy)
  
  # * optimization ------------------------------------------------------------
  
  if(print.level > 1){
    est.rez.left <- floor( (max.name.length+42-18) / 2 )
    est.rez.right <- max.name.length+42-18 - est.rez.left
    cat("\n",rep("-", est.rez.left)," The main model: ",rep("-", est.rez.right),"\n\n", sep ="")
  }
  
  
  time.05 <- proc.time()
  
  
  # mlmaximize ----------------------------------------------------------------------
  if (technique == "nr"){
    if (distribution == "e"){
      # . exponential -----------------------------------------------------------
      tymch <- tryCatch(
        .mlmaximize(
          theta0, ll = .ll.sn.exp, 
          gr = .gr.sn.exp, hess = .hess.sn.exp, gr.hess = .gr.hess.sn.exp,
          alternate = NULL, BHHH = FALSE, level = 0.99, 
          step.back = .Machine$double.eps^.5,
          reltol =  reltol, lmtol =  lmtol, 
          steptol =  .Machine$double.eps^2,
          digits = 4, when.backedup = sqrt(.Machine$double.eps), 
          max.backedup = 7,
          only.maximize = FALSE, maxit = maxit, myprod = myprod,
          y=Y, x=X, zsv=zsv, zsk=zsk, zsu=zsu, 
          n.ids=n, k=k, ksv=ksv, ksk=ksk, ksu=ksu,
          print.level = print.level,
          threads = threads),
        error = function(e) e )
    } else if (distribution == "h") {
      # . half-normal -------------------------------------------------------------
      tymch <- tryCatch(
        .mlmaximize(
          theta0, ll = .ll.sn.tn,
          gr = .gr.sn.tn, hess = .hess.sn.tn, gr.hess = .gr.hess.sn.tn,
          alternate = NULL, BHHH = FALSE, level = 0.99,
          step.back = .Machine$double.eps^.5,
          reltol =  reltol, lmtol =  lmtol, 
          steptol =  .Machine$double.eps^2,
          digits = 4, when.backedup = sqrt(.Machine$double.eps), 
          max.backedup = 7,
          only.maximize = FALSE, maxit = maxit, myprod = myprod,
          y=Y, x=X, zsv=zsv, zsk=zsk, zsu=zsu,
          n.ids=n, k=k, ksv=ksv, ksk=ksk, ksu=ksu,
          tn = tn,
          print.level = print.level,
          threads = threads),
        error = function(e) e )
    } else {
      # . truncated-normal -------------------------------------------------------------
      tymch <- tryCatch(
        .mlmaximize(
          theta0, ll = .ll.sn.tn,
          gr = .gr.sn.tn, hess = .hess.sn.tn, gr.hess = .gr.hess.sn.tn,
          alternate = NULL, BHHH = FALSE, level = 0.99,
          step.back = .Machine$double.eps^.5,
          reltol =  reltol, lmtol =  lmtol, 
          steptol =  .Machine$double.eps^2,
          digits = 4, when.backedup = sqrt(.Machine$double.eps), 
          max.backedup = 7,
          only.maximize = FALSE, maxit = maxit, myprod = myprod,
          y=Y, x=X, zsv=zsv, zsk=zsk, zsu=zsu, zmu=zmu,
          n.ids=n, k=k, ksv=ksv, ksk=ksk, ksu=ksu, kmu=kmu,
          tn = tn,
          print.level = print.level,
          threads = threads),
        error = function(e) e )
    }
  } else if (technique == "bhhh") {
    # . bhhh -------------------------------------------------------------
    
    if (distribution == "e"){
      cat("I am here\n\n")
      tymch <- tryCatch(
        .mlmaximize(
          theta0, ll = .ll.sn.exp, 
          gr = .gr.sn.exp, hess = .hess.sn.exp.bhhh, gr.hess = .gr.hess.sn.exp.bhhh,
          alternate = NULL, BHHH = FALSE, level = 0.99, 
          step.back = .Machine$double.eps^.5,
          reltol =  reltol, lmtol =  lmtol, steptol =  .Machine$double.eps^2,
          digits = 4, when.backedup = sqrt(.Machine$double.eps), max.backedup = 7,
          only.maximize = FALSE, maxit = maxit, myprod = myprod,
          y=Y, x=X, zsv=zsv, zsk=zsk, zsu=zsu, 
          n.ids=n, k=k, ksv=ksv, ksk=ksk, ksu=ksu,
          print.level = print.level,
          threads = threads),
        error = function(e) e )
    } else if (distribution == "t"){
      t
    } else {
      h
    }
  } else if (technique == "Nelder-Mead" | technique == "BFGS" | technique == "CG") {
    # optim -------------------------------------------------------------------
    # cat.print(theta0)
    if (distribution == "e"){
      # . exponential -----------------------------------------------------------
      tymch <- tryCatch(
        optim(
          theta0,
          fn = .ll.sn.exp, gr = .gr.sn.exp, 
          myprod = myprod,
          method = technique, hessian = FALSE,
          control = list(fnscale = -1, 
                         trace = optim.trace, 
                         REPORT = optim.report, 
                         maxit = ifelse(technique == 'Nelder-Mead', max(2e5, maxit), maxit), 
                         reltol = reltol),
          y=Y, x=X, zsv=zsv, zsk=zsk, zsu=zsu, 
          n.ids=n, k=k, ksv=ksv, ksk=ksk, ksu=ksu,
          threads = threads),
        error = function(e) e )
      # cat("1\n")
      if(inherits(tymch, "error")){
        if(print.level > 2){
          print(tymch)
        }
        stop("'optim' did not converge")
      }
      # print(tymch)
      if(vcetype == 'aim'){
        tymch$hessian <-
          .hess.sn.exp(tymch$par, myprod = myprod,
                       y=Y, x=X, zsv=zsv, zsk=zsk, zsu=zsu,
                       n.ids=n, k=k, ksv=ksv, ksk=ksk, ksu=ksu,
                       threads = threads)
      }
      if(vcetype == 'opg'){
        tymch$hessian <-
          .hess.sn.exp.bhhh(tymch$par, myprod = myprod,
                            y=Y, x=X, zsv=zsv, zsk=zsk, zsu=zsu,
                            n.ids=n, k=k, ksv=ksv, ksk=ksk, ksu=ksu,
                            threads = threads)
      }
      tymch$gr <- .gr.sn.exp(tymch$par, myprod = myprod,
                             y=Y, x=X, zsv=zsv, zsk=zsk, zsu=zsu,
                             n.ids=n, k=k, ksv=ksv, ksk=ksk, ksu=ksu,
                             threads = threads)
      # tymch$gr_ <- .gr.sn.exp.by.i(tymch$par, myprod = myprod,
      #                        y=Y, x=X, zsv=zsv, zsk=zsk, zsu=zsu,
      #                        n.ids=n, k=k, ksv=ksv, ksk=ksk, ksu=ksu,
      #                        threads = threads)
      # cat.print(tymch$gr)
      # cat.print(colSums(tymch$gr_))
      # cat.print(tymch$hessian)
      # tymch$vcov <- tryCatch(solve(-tymch$hessian), tol = .Machine$double.xmin * 10, error = function(e) e )
      # cat.print(tymch$vcov)
      # # ridge if needed 0
      # if(inherits(tymch$vcov, "error")){
      #   cat('I am in \n starting ridging \n')
      #   ridge <- rep(0, length(tymch$par))
      #   epsil <- 1e-6
      #   tymch$hessian1 <- tymch$hessian
      #   while(inherits(tymch$vcov, "error")){
      #     ridge <- ridge + epsil
      #     tymch$hessian1 <- tymch$hessian1 + diag(ridge, length(tymch$par))
      #     tymch$vcov <- tryCatch(solve(-tymch$hessian1), tol = .Machine$double.xmin * 10, error = function(e) e )
      #     cat.print(ridge)
      #     cat.print(inherits(tymch$vcov, "error"))
      #   }
      #   cat.print(tymch$vcov)
      # }
      # # ridge if needed 1
    } else if (distribution == "h"){
      #   . half-normal -------------------------------------------------------------
      tymch <- tryCatch(
        optim(
          theta0,
          fn = .ll.sn.tn, gr = .gr.sn.tn, hess = .hess.sn.tn,
          myprod = myprod,
          method = technique, hessian = FALSE,
          control = list(fnscale = -1,
                         trace = optim.trace,
                         REPORT = optim.report,
                         maxit = ifelse(technique == 'Nelder-Mead', max(2e5, maxit), maxit),
                         reltol = reltol),
          y=Y, x=X, zsv=zsv, zsk=zsk, zsu=zsu,
          n.ids=n, k=k, ksv=ksv, ksk=ksk, ksu=ksu,
          tn = tn,
          threads = threads),
        error = function(e) e )
      # cat("2\n")
      # print(tymch)
      tymch$hessian <- 
        .hess.sn.tn(tymch$par, myprod = myprod,
                    y=Y, x=X, zsv=zsv, zsk=zsk, zsu=zsu,
                    n.ids=n, k=k, ksv=ksv, ksk=ksk, ksu=ksu,
                    tn = tn,
                    threads = threads)
    } else  {
      # . truncated-normal -------------------------------------------------------------
      tymch <- tryCatch(
        optim(
          theta0,
          fn = .ll.sn.tn, gr = .gr.sn.tn, hess = .hess.sn.tn,
          myprod = myprod,
          method = technique, hessian = FALSE,
          control = list(fnscale = -1,
                         trace = optim.trace,
                         REPORT = optim.report,
                         maxit = ifelse(technique == 'Nelder-Mead', max(2e5, maxit), maxit),
                         reltol = reltol),
          y=Y, x=X, zsv=zsv, zsk=zsk, zsu=zsu, zmu=zmu,
          n.ids=n, k=k, ksv=ksv, ksk=ksk, ksu=ksu, kmu=kmu,
          tn = tn,
          threads = threads),
        error = function(e) e )
      # cat("3\n")
      # print(tymch)
      tymch$hessian <- 
        .hess.sn.tn(tymch$par, myprod = myprod,
                    y=Y, x=X, zsv=zsv, zsk=zsk, zsu=zsu, zmu=zmu,
                    n.ids=n, k=k, ksv=ksv, ksk=ksk, ksu=ksu, kmu=kmu,
                    tn = tn,
                    threads = threads)
    }
    
    # end of all that are done by `optim`
    if(!inherits(tymch, "error")){
      # . VCOV -----------------------------------------------------------------
      tymch$vcov <- tryCatch(solve(-tymch$hessian), tol = .Machine$double.xmin * 10, error = function(e) e )
      # cat.print(tymch$vcov)
      #  ..ridge if needed  -----------------------------------------------------
      # cat.print(tymch$par)
      # cat.print(tymch$par)
      if(inherits(tymch$vcov, "error")){
        ridge <- rep(0, length(tymch$par))
        epsil <- sqrt(.Machine$double.eps)  #1e-8
        if(print.level > 2){
          cat("\n Can't invert the hessian; starting the ridging with epsilon =",epsil,"\n")
        }
        tymch$hessian1 <- tymch$hessian
        while(inherits(tymch$vcov, "error")){
          ridge <- ridge + epsil
          tymch$hessian1 <- tymch$hessian1 + diag(ridge, length(tymch$par))
          # tymch$vcov <- tryCatch(solve(-tymch$hessian1), tol = .Machine$double.xmin * 10, error = function(e) e )
          tymch$vcov <- tryCatch(solve(-tymch$hessian1), error = function(e) e )
          # cat.print(ridge)
          # cat.print(inherits(tymch$vcov, "error"))
          if(ridge[1] > epsil*499) break
        }
        if(print.level > 2){
          cat(' The ridge that helps invery the hessian is',ridge[1],'\n')
        }
        # cat.print(tymch$vcov)
      }
    } else {
      stop("'optim' did not converge (this is repetition)")
    }
    # ridge if needed 1
  } #else {
  #stop("optimization technique is invalid")
  #}
  time.06 <- proc.time()
  est.time.sec <- (time.06-time.05)[3]
  names(est.time.sec) <- "sec"
  if(print.level >= 2){
    .timing(est.time.sec, "\nLog likelihood maximization completed in\n")
    # cat("___________________________________________________\n")
  }
  
  # time.05 <- proc.time()
  # if(optim){
  #   if(distribution == "e"){
  #     tymch <- tryCatch(
  #       optim(
  #         theta0,
  #         fn = .ll.sn.exp, gr = .gr.sn.exp, hess = .hess.sn.exp, 
  #         myprod = myprod,
  #         method = optim.method, hessian = TRUE,
  #         control = list(fnscale = -1, 
  #                        trace = trace, 
  #                        REPORT = report, 
  #                        maxit = maxit, 
  #                        reltol = reltol),
  #         y=Y, x=X, zsv=zsv, zsk=zsk, zsu=zsu, 
  #         n.ids=n, k=k, ksv=ksv, ksk=ksk, ksu=ksu,
  #         threads = threads),
  #       error = function(e) e )
  #   } else if (distribution == "t"){
  #     tymch <- tryCatch(
  #       optim(
  #         theta0,
  #         fn = .ll.sn.tn, gr = .gr.sn.tn, hess = .hess.sn.tn, 
  #         myprod = myprod,
  #         method = optim.method, hessian = TRUE,
  #         control = list(fnscale = -1, 
  #                        trace = trace, 
  #                        REPORT = report, 
  #                        maxit = maxit, 
  #                        reltol = reltol),
  #         y=Y, x=X, zsv=zsv, zsk=zsk, zsu=zsu, zmu=zmu,
  #         n.ids=n, k=k, ksv=ksv, ksk=ksk, ksu=ksu, kmu=kmu,
  #         threads = threads),
  #       error = function(e) e )
  #   } else {
  #     tymch <- tryCatch(
  #       optim(
  #         theta0,
  #         fn = .ll.sn.hn, gr = .gr.sn.hn, hess = .hess.sn.hn, 
  #         myprod = myprod,
  #         method = optim.method, hessian = TRUE,
  #         control = list(fnscale = -1, 
  #                        trace = trace, 
  #                        REPORT = report, 
  #                        maxit = maxit, 
  #                        reltol = reltol),
  #         y=Y, x=X, zsv=zsv, zsk=zsk, zsu=zsu, 
  #         n.ids=n, k=k, ksv=ksv, ksk=ksk, ksu=ksu,
  #         threads = threads),
  #       error = function(e) e )
  #   }
  # } else {
  #   if(distribution == "e"){
  #     tymch <- tryCatch(
  #       .mlmaximize(
  #         theta0, ll = .ll.sn.exp, 
  #         gr = .gr.sn.exp, hess = .hess.sn.exp, gr.hess = .gr.hess.sn.exp,
  #         alternate = NULL, BHHH = FALSE, level = 0.99, 
  #         step.back = .Machine$double.eps^.5,
  #         reltol =  reltol, lmtol =  lmtol, steptol =  .Machine$double.eps^2,
  #         digits = 4, when.backedup = sqrt(.Machine$double.eps), max.backedup = 7,
  #         only.maximize = FALSE, maxit = maxit, myprod = myprod,
  #         y=Y, x=X, zsv=zsv, zsk=zsk, zsu=zsu, 
  #         n.ids=n, k=k, ksv=ksv, ksk=ksk, ksu=ksu,
  #         print.level = print.level,
  #         threads = threads),
  #       error = function(e) e )
  #   } else if (distribution == "t"){
  #     tymch <- tryCatch(
  #       .mlmaximize(
  #         theta0, ll = .ll.sn.tn, 
  #         gr = .gr.sn.tn, hess = .hess.sn.tn, gr.hess = .gr.hess.sn.tn,
  #         alternate = NULL, BHHH = FALSE, level = 0.99, 
  #         step.back = .Machine$double.eps^.5,
  #         reltol =  reltol, lmtol =  lmtol, steptol =  .Machine$double.eps^2,
  #         digits = 4, when.backedup = sqrt(.Machine$double.eps), max.backedup = 7,
  #         only.maximize = FALSE, maxit = maxit, myprod = myprod,
  #         y=Y, x=X, zsv=zsv, zsk=zsk, zsu=zsu, zmu=zmu,
  #         n.ids=n, k=k, ksv=ksv, ksk=ksk, ksu=ksu, kmu=kmu,
  #         print.level = print.level,
  #         threads = threads),
  #       error = function(e) e )
  #   } else {
  #     tymch <- tryCatch(
  #       .mlmaximize(
  #         theta0, ll = .ll.sn.exp, 
  #         gr = .gr.sn.hn, hess = .hess.sn.hn, gr.hess = .gr.hess.sn.hn,
  #         alternate = NULL, BHHH = FALSE, level = 0.99, 
  #         step.back = .Machine$double.eps^.5,
  #         reltol =  reltol, lmtol =  lmtol, steptol =  .Machine$double.eps^2,
  #         digits = 4, when.backedup = sqrt(.Machine$double.eps), max.backedup = 7,
  #         only.maximize = FALSE, maxit = maxit, myprod = myprod,
  #         y=Y, x=X, zsv=zsv, zsk=zsk, zsu=zsu, 
  #         n.ids=n, k=k, ksv=ksv, ksk=ksk, ksu=ksu,
  #         print.level = print.level,
  #         threads = threads),
  #       error = function(e) e )
  #   }
  # }
  # time.06 <- proc.time()
  # est.time.sec <- (time.06-time.05)[3]
  # names(est.time.sec) <- "sec"
  # if(print.level >= 2){
  #   .timing(est.time.sec, "\nLog likelihood maximization completed in\n")
  #   # cat("___________________________________________________\n")
  # }
  
  # print(m1$table)
  
  if(inherits(tymch, "error")){
    if(print.level > 0){
      print(tymch)
      cat(" There was an error in optimizer \n")
    }
  }
  # if(optim){
  #   tymch$vcov <- tryCatch(solve(-tymch$hessian), error = function(e) e )
  # }
  
  if (technique == "Nelder-Mead" | technique == "BFGS" | technique == "CG") {
    tymch$ll <- tymch$value
    tymch$value <- NULL
  }
  
  # cat.print(inherits(tymch$vcov, "error"))
  if(inherits(tymch$vcov, "error")){
    if(print.level > 2){
      cat("don't know what to do\n")
    }
    class(tymch) <- "snreg"
    # return(tymch)
  } else {
    if (technique == "Nelder-Mead" | technique == "BFGS" | technique == "CG"){
      if(print.level > 2){
        cat(paste("\nFinal log likelihood = ",format(tymch$ll, digits = 13),"\n", sep = ""), sep = "")
      }
    }
    tymch$sds <- suppressWarnings( sqrt(diag(tymch$vcov)) )
    tymch$ctab <- cbind(tymch$par, tymch$sds, tymch$par/tymch$sds,2*(1-pnorm(abs(-tymch$par/tymch$sds))))
    colnames(tymch$ctab) <- c("Estimate", "Std.Err", "Z value", "Pr(>z)")
    colnames(tymch$ctab) <- c("Coef.", "SE ", "z ",  "P>|z|")
    coef.mat.print(tymch$ctab)
    tymch$ctab2d2 <- coef.mat.print(tymch$ctab, digits = 2)
    tymch$ctab2d3 <- coef.mat.print(tymch$ctab, digits = 3)
    (tymch$ctab2d4 <- coef.mat.print(tymch$ctab, digits = 4) )
    # cat("\n")
    # print(tymch$ctab2d4)
  }
  
  # output ------------------------------------------------------------------
  
  output <- tymch$ctab#cbind(round(par, digits = digits), round(se, digits = digits), round(par/se,digits = 2), round(pnorm(abs(par/se), lower.tail = FALSE)*2, digits = digits))
  colnames(output) <- c("Coef.", "SE ", "z ",  "P>|z|")
  
  max.name.length <- max(nchar(names(theta0) ))
  
  if(print.level >= 1){
    cat("",rep("_", max.name.length+42-1),"", "\n", sep = "")
    if(prod){
      cat("\nCross-sectional stochastic (production) frontier model\n",  sep = "")}
    else {
      cat("\nCross-sectional stochastic (cost) frontier model\n",  sep = "")
    }
    cat("\nDistributional assumptions\n\n", sep = "")
    Assumptions <- rep("heteroskedastic",2)
    if(ksv==1){
      Assumptions[1] <- "homoskedastic"
    }
    if(ksu==1){
      Assumptions[2] <- "homoskedastic"
    }
    Distribution = c("skew normal ", "half-normal ")
    if(distribution == "t"){
      Distribution[2] <- "truncated-normal "
    }   
    if(distribution == "e"){
      Distribution[2] <- "exponential "
    }
    a1 <- data.frame(
      Component = c("Random noise: ","Inefficiency: "),
      Distribution = Distribution,
      Assumption = Assumptions
    )
    print(a1, quote = FALSE, right = FALSE)
    cat("\nNumber of observations = ", n, "\n", sep = "")
    # max.name.length <- max(nchar(row.names(output)))
    est.rez.left <- floor( (max.name.length+42-22) / 2 )
    est.rez.right <- max.name.length+42-22 - est.rez.left
    cat("\n",rep("-", est.rez.left)," Estimation results: ",rep("-", est.rez.right),"\n\n", sep ="")
    # cat("\n--------------- Estimation results: --------------\n\n", sep = "")
    .printoutcs(output, digits = digits, k = k, ksv = ksv, ksk = ksk, ksu = ksu, kmu = kmu, na.print = "NA", dist = distribution, max.name.length = max.name.length)
  }
  
  # cat('Efficiency \n')
  # efficiency --------------------------------------------------------------
  
  if(distribution == "e"){
    u.tymch <- .u.sn.exp(
      theta = tymch$par, myprod = myprod,
      y=Y, x=X, zsv=zsv, zsk=zsk, zsu=zsu,
      n.ids=n, k=k, ksv=ksv, ksk=ksk, ksu=ksu,
      threads = threads)
    
    u     <- u.tymch$u
    u.tymch$eps -> eps0
    sv <- u.tymch$sv
    al0 <- u.tymch$al
    lam <- u.tymch$lam
    
    # eps0  <- Y - X %*%    tymch$par[1                        :k]
    # cat.print( summary(u.tymch$eps - eps0) )
    # sv    <- sqrt(exp( zsv %*% tymch$par[(k + 1)                  :(k + ksv)]))
    # cat.print(summary( sv - u.tymch$sv))
    # al0   <- zsk %*%      tymch$par[(k + ksv + 1)            :(k + ksv + ksk)]
    # cat.print(summary( al0 - u.tymch$al))
    # lam   <- 1 / sqrt(exp( zsu %*% tymch$par[(k + ksv + ksk  + 1)     :(k + ksv + ksk + ksu)] ))
    # cat.print(summary( lam - u.tymch$lam))
    
    
    # epsr  <- eps0 + sv * sqrt(2/pi) * al0 / sqrt(1 + al0^2)
    # u1    <- (myprod*epsr + lam*sv^2)/sv
    # b1    <- al0 * myprod
    # b2    <- sqrt(1+b1^2)
    # a1    <- -b1 * lam * sv
    # a2    <- a1 / b2
    # 
    # term1 <- -TOwen( u1, a2 / u1 )
    # term2 <- -TOwen( a2, u1 / a2 )
    # term3 <-  TOwen( u1, b1 + a1 / u1 )
    # term4 <-  TOwen( a2, b1 + u1 * (1 + b1^2) / a1 )
    # term5 <- pnorm(a2) * pnorm(-u1)
    # t122  <- term1 + term2 + term3 + term4 + term5
    # 
    # t135  <- b1/b2*dnorm(a2)*pnorm(-u1*b2 - b1*a2) +
    #   dnorm(u1)*pnorm(a1+b1*u1)
    # cat.print(summary(t135))
    # cat.print(summary(t122))
    # cat.print(summary(u1))
    # cat.print(summary(sv))
    # cat.print(summary(t135 / t122))
    # 
    # u_     <- sv * ( t135 / t122 - u1)
    # 
    # # cat.print(summary(u))
    # # cat.print(summary(u_-u))
    # tymch2 <- data.frame(round(eps0,3) , round(epsr,3), round(sv,3), round(al0,3), round(lam,3) , round(a1,3), round(b1,3), round(u1,3), term1,term2,term3,term4,term5,t122, t135, u_, u,term2+term4,u1 / a2, b1 + u1 * (1 + b1^2) / a1 )
    # colnames(tymch2) <- c("eps0","epsr","sv","al0","lam","a1","b1","u1","term1","term2","term3","term4","term5","t122","t135","u_","u","term2-4","u1-a2","b1ldots")
    # tymch2 <- tymch2[order(tymch2$u_),]
    # cat.print(tymch2)
    # cat.print(tymch2[u_<0,])
  } else {
    eps0  <- Y - X %*%    tymch$par[1                        :k]
    sv2   <- exp( zsv %*% tymch$par[(k + 1)                  :(k + ksv)])
    al0   <- zsk %*%      tymch$par[(k + ksv + 1)            :(k + ksv + ksk)]
    su2   <- exp( zsu %*% tymch$par[(k + ksv + ksk  + 1)     :(k + ksv + ksk + ksu)] )
    if(distribution == "h"){
      mu0 <- 0
    } else {
      mu0   <- zmu %*%      tymch$par[(k + ksv + ksk + ksu + 1):(k + ksv + ksk + ksu + kmu)]
    }
    
    sv    <- sqrt(sv2)
    su    <- sqrt(su2)
    sig   <- sqrt(sv2 + su2)
    sstar <- sv*su/sig
    epsr  <- eps0 + sv * sqrt(2/pi) * al0 / sqrt(1+al0^2)
    mu1   <- (mu0*sv2 - epsr*myprod*su2) / sig^2
    b1    <- al0 / sv * myprod * sstar 
    b2    <- sqrt(1+b1^2)
    a1    <- al0 / sv * (epsr + myprod * mu1)
    a2    <- a1 / b2
    u1    <- -mu1 / sstar
    
    term1 <- -TOwen( u1, a2 / u1 )
    term2 <- -TOwen( a2, u1 / a2 )
    term3 <-  TOwen( u1, b1 + a1 / u1 )
    term4 <-  TOwen( a2, b1 + u1 * (1 + b1^2) / a1 )
    term5 <- pnorm(a2) * pnorm(-u1)
    t92   <- term1 + term2 + term3 + term4 + term5
    
    t100  <- b1/b2*dnorm(a1/b2)*pnorm(-u1*b2 - a1*b1/b2) +
      dnorm(u1)*pnorm(a1+b1*u1)
    
    u     <- mu1 + sstar * t100 / t92
  }
  # cat('Auxiliary parameters \n')
  # Auxiliary parameters
  
  tymch $ resid   <- as.vector(unname(eps0))
  shat2           <- var( tymch$resid )
  tymch $ shat2   <- shat2
  tymch $ shat    <- sqrt(shat2)
  tymch $ RSS     <- crossprod(eps0)
  tymch $ aic     <- log((n-1)/n*shat2)+1+2*(k.all+1)/n
  tymch $ bic     <- log((n-1)/n*shat2)+1+(k.all+1)*log(n)/n
  tymch $ aic     <- 2*k.all - 2*tymch$ll
  tymch $ bic     <- log(n)*k.all - 2*tymch$ll
  tymch $ Mallows <- tymch $ RSS/shat2 - n + 2*k.all
  tymch $ coef    <- tymch$par
  tymch $ esample <- esample
  tymch $ sv      <- as.vector(unname(sv))
  tymch $ n       <- n
  if(distribution == "e"){
    tymch$ su    <- 1/as.vector(unname(lam))
  } else {
    tymch$ su    <- unname(su)
  }
  tymch$ skewness   <- as.vector(unname(al0))
  tymch$ u     <- as.vector(unname(u))
  tymch$ eff   <- exp(-tymch$ u)
  
  # cat.print(tymch$coef)
  # cat.print(tymch$par)
  # cat.print(tymch$vcov)
  # cat.print(coef.names.full)
  
  names(tymch$coef) <- names(tymch$par) <- colnames(tymch$vcov) <- rownames(tymch$vcov) <-
    coef.names.full
  
  # tymch2 <- data.frame(t135, t122, t135 / t122, u1, u, tymch$ eff)
  # colnames(tymch2) <- c("t135", " t122", " t135 / t122", " e_rs", " u", " te")
  # cat.print(tymch2[tymch$ eff > .999 | tymch$ eff < .3,])
  # cat.print(summary(tymch$eff))
  # cat("\n")
  # cat.print(summary(tymch2))
  
  # if(print.level >= 3) cat.print(summary(tymch$ eff))
  myeff <- ifelse(prod, "technical", "cost")
  if(print.level >= 3){
    eff.name <- paste0("Summary of ",myeff," efficiencies")
    len.eff.name <- nchar(eff.name)
    est.eff.left <- floor( (max.name.length+42-len.eff.name-4) / 2 )
    est.eff.right <- max.name.length+42-len.eff.name-4 - est.eff.left
    cat("\n",rep("-", est.eff.left)," ",eff.name,": ",rep("-", est.eff.right),"\n\n", sep ="")
    # eff1 <- data.frame(TE_JLMS = tymch$ eff)#eff[,2:4]
    # colnames(eff1) <- c("TE_JLMS")
    # cat.print(eff1)
    # colnames(eff1) <- formatC(colnames(eff1), width = 4, flag = "-")
    cat("",rep(" ", est.eff.left+1),"JLMS:= exp(-E[ui|ei])\n", sep = "")
    # cat("",rep(" ", est.eff.left+1),"Mode:= exp(-M[ui|ei])\n", sep = "")
    # cat("",rep(" ", est.eff.left+3),"BC:= E[exp(-ui)|ei]\n", sep = "")
    cat("\n")
    .su(tymch$ eff, transpose = TRUE, print = TRUE, names = "TE_JLMS")
    
    # cat("\n=================")
    # cat(" Summary of ",myeff," efficiencies, exp(-E[ui|ei]): \n\n", sep = "")
    # .su1(eff1[,1, drop = FALSE], transpose = TRUE)
    #
    # cat("\n=================")
    # cat(" Summary of ",myeff," efficiencies, exp(-M[ui|ei]): \n\n", sep = "")
    # .su1(eff1[,2, drop = FALSE])
    #
    # cat("\n=================")
    # cat(" Summary of ",myeff," efficiencies, E[exp(-ui)|ei]: \n\n", sep = "")
    # .su1(eff1[,3, drop = FALSE])
    # cat("\n\n")
  }
  
  class(tymch) <- "snreg"
  return(tymch)
  
}




#


