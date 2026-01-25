
#' Linear Regression with Skew-Normal Errors
#'
#' @title Linear Regression with Skew-Normal Errors
#'
#' @description
#' \code{snreg} fits a linear regression model where the disturbance term follows
#' a skew-normal distribution. The function supports multiplicative
#' heteroskedasticity of the noise variance via a log-linear specification
#' (\code{ln.var.v}) and allows the skewness parameter to vary linearly with
#' exogenous variables (\code{skew.v}).
#'
#' @param formula
#' an object of class \code{formula} specifying the regression:
#' typically \code{y ~ x1 + ...}, where \code{y} is the dependent variable
#' and \code{x}'s are regressors.
#'
#' @param data
#' an optional \code{data.frame} containing the variables in \code{formula}.
#' If not found in \code{data}, variables are taken from \code{environment(formula)}.
#'
#' @param subset
#' an optional logical or numeric vector specifying the subset of observations
#' to be used in estimation.
#'
#' @param start.sk
#' numeric. Initial value for the (global) skewness parameter of the noise;
#' can be \code{NULL} if \code{skew.v} is supplied with its own coefficients to initialize.
#'
#' @param ln.var.v
#' optional one-sided formula; e.g. \code{ln.var.v ~ z1 + z2}. Specifies
#' exogenous variables entering the (log) variance of the random noise component.
#' If \code{NULL}, the noise variance is homoskedastic.
#'
#' @param skew.v
#' optional one-sided formula; e.g. \code{skew.v ~ z3 + z4}. Specifies exogenous
#' variables determining the skewness of the noise via a linear index; if
#' \code{NULL}, the skewness is constant (scalar).
#'
#' @param start.val
#' optional numeric vector of starting values for all free parameters
#' (regression coefficients, variance/heteroskedasticity parameters, skewness parameters).
#'
#' @param technique
#' character vector giving the preferred maximization routine(s) in order of
#' preference. Currently recognized keywords include \code{"nr"} (Newton–Raphson),
#' \code{"bhhh"}, \code{"nm"} (Nelder–Mead), \code{"bfgs"}, \code{"cg"}.
#' This scaffold does not implement them yet, but records the choice.
#'
#' @param vcetype
#' character specifying the variance-covariance estimator type:
#' \code{"aim"} for the approximated information matrix or \code{"opg"}
#' for the outer product of gradients. Default is \code{"aim"}.
#'
#' @param lmtol
#' numeric. Convergence tolerance based on the scaled gradient (if applicable).
#' Default is \code{1e-5}.
#'
#' @param reltol
#' numeric. Relative convergence tolerance for likelihood maximization.
#' Default is \code{1e-12}.
#'
#' @param maxit
#' integer. Maximum number of iterations for the optimizer. Default is \code{199}.
#'
#' @param report
#' integer. Verbosity for reporting progress (if implemented). Default is \code{1}.
#'
#' @param trace
#' integer. If positive, tracing information is printed (if implemented).
#' Default is \code{1}.
#'
#' @param print.level
#' integer. Printing level for summaries: \code{1}—print estimation results;
#' \code{2}—print optimization details; \code{3}—print compact summary. Default \code{3}.
#'
#' @param digits
#' integer. Number of digits for printing. Default \code{4}.
#'
#' @param threads
#' integer. Number of threads (placeholder for parallel implementations).
#' Default \code{1}.
#'
#' @param only.data
#' logical. If \code{TRUE}, the function returns only the constructed model
#' matrices and design sets (no estimation). Default \code{FALSE}.
#'
#' @param ...
#' additional arguments reserved for future methods (e.g., box constraints).
#'

#' @details
#' The model is
#' \deqn{y_i = x_i^\top \beta + \varepsilon_i,\quad \varepsilon_i \sim SN(0, \sigma_i^2, \alpha_i),}
#' where \eqn{SN} denotes the skew-normal distribution (Azzalini).
#'
#' Heteroskedasticity in the noise variance (if specified via \code{ln.var.v}) is modeled as
#' \deqn{\log(\sigma_i^2) = w_i^\top \gamma_v,}
#' and the (optional) covariate-driven skewness (if specified via \code{skew.v}) as
#' \deqn{\alpha_i = s_i^\top \delta.}
#'
#' This function constructs the model frame and design matrices for
#' \eqn{\beta}, \eqn{\gamma_v}, and \eqn{\delta}, and is designed to be paired with
#' a maximum likelihood routine to estimate parameters and (optionally) their
#' asymptotic covariance via either AIM or OPG.
#'
#' @return
#' An object of class \code{"snreg"} with elements:
#' \itemize{
#'   \item{\code{call}}{ — matched function call.}
#'   \item{\code{terms}}{ — model terms for the main regression.}
#'   \item{\code{model}}{ — list with constructed data: \code{y}, \code{X}, \code{Zv}, \code{Zs}.}
#'   \item{\code{coef}}{ — named vector of MLEs (placeholder \code{numeric(0)} in scaffold).}
#'   \item{\code{vcov}}{ — variance–covariance matrix (placeholder).}
#'   \item{\code{loglik}}{ — log-likelihood at the solution (placeholder \code{NA}).}
#'   \item{\code{esample}}{ — logical vector indicating the estimation sample.}
#'   \item{\code{controls}}{ — list of control parameters and settings.}
#' }
#'
#' @references
#' Azzalini, A. (1985).
#' \emph{A Class of Distributions Which Includes the Normal Ones}.
#' Scandinavian Journal of Statistics, 12(2), 171–178.
#'
#' Azzalini, A., & Capitanio, A. (2014).
#' \emph{The Skew-Normal and Related Families}.
#' Cambridge University Press.
#'
#' @examples
#' \dontrun{
#'   # Simulated usage (replace with real data)
#'   set.seed(1)
#'   n <- 200
#'   x1 <- rnorm(n); x2 <- runif(n)
#'   X  <- cbind(1, x1, x2)
#'   beta <- c(1, 2, -1)
#'   y <- X %*% beta + rnorm(n)
#'   df <- data.frame(y = as.numeric(y), x1 = x1, x2 = x2,
#'                    z1 = rnorm(n), z2 = rnorm(n))
#'
#'   # Constant variance, constant skewness
#'   m0 <- snreg(
#'     formula = y ~ x1 + x2, data = df,
#'     ln.var.v = NULL, skew.v = NULL,
#'     technique = c("nr")
#'   )
#'   str(m0)
#'
#'   # Heteroskedastic variance and variable skewness
#'   m1 <- snreg(
#'     formula = y ~ x1 + x2, data = df,
#'     ln.var.v = ~ z1 + z2,
#'     skew.v   = ~ z1,
#'     technique = c("bfgs"), vcetype = "opg"
#'   )
#'   str(m1)
#' }
#'
#' @keywords regression skew-normal heteroskedasticity maximum-likelihood
#' @export
snreg <- function (
  formula, data, subset,
  start.sk  = NULL,#.5*(2*prod-1),
  ln.var.v  = NULL,
  skew.v    = NULL,
  start.val = NULL,
  technique = c('nr'),   #,'bhhh','nm', 'bfgs', 'cg'),
  vcetype   = c('aim'),  #, 'opg'), # `approximated information matrix` or `outer product of gradients`
  lmtol     = 1e-5,  
  reltol    = 1e-12, 
  maxit     = 199, 
  report    = 1,
  trace     = 1,
  print.level = 3, 
  digits    = 4,
  threads   = 1,
  only.data = FALSE,
  ...){
  
  # threads <- 1
  
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

  form1 <- as.Formula(formula, ln.var.v, skew.v)

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

  # from `sf`
  
  colnames(X)[1]   <- "(Intercept)"
  colnames(zsv)[1] <- "(Intercept)"
  colnames(zsk)[1] <- "(Intercept)"

  abbr.length <- 34
  names_x <- abbreviate(colnames(X), abbr.length+6, strict = TRUE, dot = FALSE, method = "both.sides")
  # names_zu <-  abbreviate(colnames(Zu), 9, strict = TRUE, dot = FALSE)
  # names_zv <-  abbreviate(colnames(Zv), 9, strict = TRUE, dot = FALSE)
  # Zu_colnames <- abbreviate(colnames(Zu), abbr.length, strict = TRUE, dot = FALSE, method = "both.sides")
  names_zv <- paste0("lnVARv0i_", abbreviate(colnames(zsv), abbr.length, strict = TRUE, dot = FALSE, method = "both.sides"))
  names_zk <- paste0("Skew_v0i_", abbreviate(colnames(zsk), abbr.length, strict = TRUE, dot = FALSE, method = "both.sides"))
  

  coef.names.full <- c(
    names_x,
    paste("lnVARv0i_",c("(Intercept)", names_zv[-1]),"", sep = ""),
    paste("Skew_v0i_",c("(Intercept)", names_zk[-1]),"", sep = "")
  )
  # cat.print(coef.names.full)
  # coef.names.full <- c(
  #   names_x,  names_zv, names_zk
  # )

  # starting values ---------------------------------------------------------
  
  # OLS
  tymch2 <- lm(Y ~ X - 1)
  olsres  <- resid(tymch2)
  shat <- mean(olsres^2)
  
  gamma1.0 <- min(.99, mean(olsres^3)/shat )
  
  r  <- sign(gamma1.0)*(2*abs(gamma1.0)/(4-pi))^(1/3)
  delta <- r/(sqrt(2/pi) *sqrt(1+r^2))
  if(abs(delta) > 1){
    alpha.0 <- delta/sqrt(1-delta^2)
  } else {
    alpha.0 <- sign(delta) * 3
  }
  
  # cat.print(alpha.0)
  
  # b <- sqrt(2/pi) 
  # sd <- 1.3069693
  # gamma1 <- 0.6670236
  # r  <- sign(gamma1)*(2*abs(gamma1)/(4-pi))^(1/3)
  # delta <- r/(b*sqrt(1+r^2))
  # ( alpha <- delta/sqrt(1-delta^2) )
  # (mu.z <- b*delta )
  # (sd.z <- sqrt(1-mu.z^2) )
  # ( omega <- sd/sd.z )

  if(is.null(start.sk)){
    # sk.initial <- tymch1$par[(k+ksv+1):(k+ksv+ksk)] * .5
    if(ksk==1){
      sk.initial <- alpha.0
    } else {
      sk.initial <- c(alpha.0, rep(0,ksk-1))
    }
  } else {
    if(ksk==1){
      sk.initial <- start.sk
    } else {
      sk.initial <- c(start.sk, rep(0,ksk-1))
    }
  }
  


  lnsv00 <- log(mean(olsres^2))
  if(ksv==1){
    sv.initial <- lnsv00
  } else {
    sv.initial <- c(lnsv00, rep(0,ksv-1))
  }
  
  theta0 <- c( coef(tymch2), sv.initial, sk.initial)
  names(theta0) <- coef.names.full
  
  if(print.level > 0) {
    cat.print(theta0)
  }
  
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
    trace1 <- 10
  }
  
  # * optimization ------------------------------------------------------------
  
  if(technique == 'BFGS'){
    tymch <- optim(
      par = theta0, fn = .ll.sn, gr = .gr.sn, y = Y, x = X, zsv = zsv, zsk = zsk, k = k, ksv = ksv, ksk = ksk,
      method = technique, control = list(fnscale = -1, trace = trace1, maxit = 10000), 
      hessian = FALSE)
    
    if(vcetype == 'opg'){
      gri <- .gr.sn.by.i(tymch$par, y = Y, x = X, zsv = zsv, zsk = zsk, k = k, ksv = ksv, ksk = ksk)
      tymch$hessian <- -crossprod(gri)
    } else {
      tymch$hessian <- hessian2(funcg = .gr.sn, at = tymch$par, y = Y, x = X, zsv = zsv, zsk = zsk, k = k, ksv = ksv, ksk = ksk)
    }
    tymch$coef <- tymch$par
    tymch$ll <- tymch$value
    tymch$vcov <- tryCatch(solve(-tymch$hessian), tol = .Machine$double.xmin * 10, error = function(e) e )
    
  } else if (technique == 'nr'){
    tymch <- .mlmaximize(theta0, ll = .ll.sn, gr = .gr.sn, hess = .hess.sn, gr.hess = .hess.gr.sn, alternate = NULL, BHHH = FALSE, level = 0.99, step.back = .Machine$double.eps^.5, reltol =  sqrt(.Machine$double.eps), lmtol =  1e-4, steptol =  .Machine$double.eps, digits = 4, when.backedup = sqrt(.Machine$double.eps), max.backedup = 7, print.level = print.level, only.maximize = FALSE, maxit = 250, 
                               y = Y, x = X, zsv = zsv, zsk = zsk, k = k, ksv = ksv, ksk = ksk                   )
    tymch$coef <- tymch$par
    # print(tymch)
  }
  
  # cat.print(names(tymch))
  
  
  
  # cat.print(tymch1)
  # table with coefs
  # tymch$par  <- c(tymch$par, sv = unname(sqrt(exp(tymch$par[k+1]))) )
  # tymch$sds  <- suppressWarnings( sqrt(diag(solve(-tymch$hessian))) )
  # tymch$sds  <- c(tymch$sds, tymch$sds[k+1]*exp(tymch$par[k+1]))
  
  tymch$sds  <- suppressWarnings( sqrt(diag(tymch$vcov )) )
  
  

  tymch$ctab <- cbind(tymch$par, tymch$sds, tymch$par/tymch$sds,2*(1-pnorm(abs(-tymch$par/tymch$sds))))
  colnames(tymch$ctab) <- c("Estimate", "Std.Err", "Z value", "Pr(>z)")
  
  if(print.level >= 1){
    printCoefmat(tymch$ctab)
    cat("",rep("_", max.name.length+42-1),"", "\n", sep = "")
  }
  # cat.print(tymch1)
  # return(tymch1)
  

  # cat('Auxiliary parameters \n')
  # Auxiliary parameters
  
  # tymch $ resid   <- as.vector(unname(eps0))
  # shat2           <- var( tymch$resid )
  # tymch $ shat2   <- shat2
  # tymch $ shat    <- sqrt(shat2)
  eps0  <- Y - X %*%    tymch$par[1                        :k]
  tymch $ RSS     <- crossprod(eps0)
  # tymch $ aic     <- log((n-1)/n*shat2)+1+2*(k.all+1)/n
  # tymch $ bic     <- log((n-1)/n*shat2)+1+(k.all+1)*log(n)/n
  # tymch $ aic     <- 2*k.all - 2*tymch$ll
  # tymch $ bic     <- log(n)*k.all - 2*tymch$ll
  # tymch $ Mallows <- tymch $ RSS/shat2 - n + 2*k.all
  # tymch $ coef    <- tymch$par
  tymch $ esample <- esample
  # tymch $ sv      <- as.vector(unname(sv))
  tymch $ n       <- n
  # if(distribution == "e"){
  #   tymch$ su    <- 1/as.vector(unname(lam))
  # } else {
  #   tymch$ su    <- unname(su)
  # }
  tymch$ skewness   <- as.vector(unname(zsk %*%      tymch$par[(k + ksv + 1)            :(k + ksv + ksk)]))
  # tymch$ u     <- as.vector(unname(u))
  # tymch$ eff   <- exp(-tymch$ u)
  
  # cat.print(tymch$coef)
  # cat.print(tymch$par)
  # cat.print(tymch$vcov)
  # cat.print(coef.names.full)
  
  # cat.print(names(tymch$coef) )
  # cat.print(names(tymch$par))
  # cat.print(colnames(tymch$vcov))
  # cat.print(rownames(tymch$vcov))
  # cat.print(coef.names.full)
  if(technique != 'nr'){
    names(tymch$coef) <- names(tymch$par) <- colnames(tymch$vcov) <- rownames(tymch$vcov) <-
      coef.names.full
  }

  
  class(tymch) <- "snreg"
  return(tymch)
  
}




#


