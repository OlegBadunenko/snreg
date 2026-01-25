
#' Linear Model by Maximum Likelihood (with optional heteroskedasticity)
#'
#' @title Linear Regression via MLE
#'
#' @description
#' \code{lm.mle} fits a linear regression model by maximum likelihood, allowing
#' for optional multiplicative heteroskedasticity in the disturbance variance via
#' a log-linear specification provided through \code{ln.var.v}.
#'
#' @param formula
#' an object of class \code{formula} specifying the regression:
#' typically \code{y ~ x1 + ...}, where \code{y} is the dependent variable and
#' the \code{x}'s are regressors.
#'
#' @param data
#' an optional \code{data.frame} containing the variables referenced in
#' \code{formula}. If not found in \code{data}, variables are taken from
#' \code{environment(formula)}.
#'
#' @param subset
#' an optional logical or numeric vector specifying the subset of observations
#' to be used in estimation.
#'
#' @param ln.var.v
#' optional one-sided formula; e.g. \code{ln.var.v ~ z1 + z2}. When provided,
#' the error variance is modeled as \eqn{\log(\sigma_i^2) = w_i^\top \gamma_v}.
#' If \code{NULL}, the variance is homoskedastic.
#'
#' @param technique
#' character vector specifying the preferred optimization routine(s) in order of
#' preference. Recognized keywords (for future implementation) include \code{"nr"}
#' (Newton–Raphson), \code{"bhhh"}, \code{"nm"} (Nelder–Mead), \code{"bfgs"},
#' and \code{"cg"}. Default is \code{"nr"}. This scaffold records but does not
#' execute the chosen routine.
#'
#' @param vcetype
#' character specifying the variance–covariance estimator type:
#' \code{"aim"} for (approximate) information matrix or \code{"opg"} for the outer
#' product of gradients. Default is \code{"aim"} (recorded; not yet computed here).
#'
#' @param lmtol
#' numeric. Convergence tolerance based on scaled gradient (when applicable).
#' Default \code{1e-5}.
#'
#' @param reltol
#' numeric. Relative convergence tolerance for likelihood maximization.
#' Default \code{1e-12}.
#'
#' @param maxit
#' integer. Maximum number of iterations for the optimizer. Default \code{199}.
#'
#' @param report
#' integer. Verbosity level for reporting progress (if implemented). Default \code{1}.
#'
#' @param trace
#' integer. Trace level for optimization (if implemented). Default \code{1}.
#'
#' @param print.level
#' integer. Printing level for summaries. Default \code{3}.
#'
#' @param digits
#' integer. Number of digits for printing. Default \code{4}.
#'
#' @param threads
#' integer. Number of threads (placeholder for parallel implementations). Default \code{8}.
#'
#' @param only.data
#' logical. If \code{TRUE}, returns only constructed data/matrices without
#' estimation. Default \code{FALSE}.
#'
#' @param ...
#' additional arguments reserved for future methods (e.g., bounds, penalties).
#'
#' @details
#' The model is
#' \deqn{y_i = x_i^\top \beta + \varepsilon_i,\quad \varepsilon_i \sim \mathcal{N}(0, \sigma_i^2).}
#' When \code{ln.var.v} is supplied, the variance follows
#' \deqn{\log(\sigma_i^2) = w_i^\top \gamma_v,}
#' otherwise \eqn{\sigma_i^2 = \sigma^2} is constant (homoskedastic).
#'
#' This scaffold:
#' \itemize{
#'   \item Builds the model frame and \code{X}, \code{y}.
#'   \item Builds \code{Zv} for the log-variance index when \code{ln.var.v} is provided.
#'   \item Returns a structured object with placeholders for \code{coef}, \code{vcov}, \code{loglik}.
#' }
#' Insert your MLE engine to estimate \eqn{\beta}, and (optionally) \eqn{\sigma^2} or
#' \eqn{\gamma_v}; compute standard errors via AIM/OPG as required by \code{vcetype}.
#'
#' @return
#' A list of class \code{"lm.mle"} containing:
#' \itemize{
#'   \item{\code{call}}{ — matched call.}
#'   \item{\code{terms}}{ — model terms.}
#'   \item{\code{model}}{ — list with constructed components: \code{y}, \code{X}, \code{Zv}.}
#'   \item{\code{coef}}{ — named numeric vector of parameter estimates (placeholder).}
#'   \item{\code{vcov}}{ — variance–covariance matrix (placeholder).}
#'   \item{\code{loglik}}{ — log-likelihood at the solution (placeholder).}
#'   \item{\code{esample}}{ — logical vector indicating the estimation sample.}
#'   \item{\code{controls}}{ — list of control parameters and settings.}
#' }
#'
#' @examples
#' \dontrun{
#'   set.seed(42)
#'   n  <- 300
#'   x1 <- rnorm(n); x2 <- runif(n)
#'   y  <- 1 + 2*x1 - 1.5*x2 + rnorm(n, sd = 0.7)
#'   df <- data.frame(y, x1, x2, z1 = rnorm(n))
#'
#'   # Homoskedastic MLE scaffold
#'   fit0 <- lm.mle(y ~ x1 + x2, data = df)
#'   str(fit0)
#'
#'   # Heteroskedastic variance via log-linear index
#'   fit1 <- lm.mle(
#'     y ~ x1 + x2, data = df,
#'     ln.var.v = ~ z1,
#'     technique = c("nr"), vcetype = "aim"
#'   )
#'   str(fit1)
#'
#'   # Design matrices only (no estimation)
#'   mats <- lm.mle(y ~ x1 + x2, data = df, ln.var.v = ~ z1, only.data = TRUE)
#'   names(mats$model)
#' }
#'
#' @keywords regression maximum-likelihood heteroskedasticity
#' @export
lm.mle <- function (
  formula, data, subset,
  ln.var.v  = NULL,
  technique = c('nr'),#,'bhhh','nm', 'bfgs', 'cg'),
  vcetype   = c('aim'),#, 'opg'), # `approximated information matrix` or `outer product of gradients`
  lmtol     = 1e-5,  
  reltol    = 1e-12, 
  maxit     = 199, 
  report    = 1,
  trace     = 1,
  print.level = 3, 
  digits    = 4,
  threads   = 8,
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

  form1 <- as.Formula(formula, ln.var.v)

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

  # from `sf`
  
  colnames(X)[1]   <- "(Intercept)"
  colnames(zsv)[1] <- "(Intercept)"

  abbr.length <- 34
  names_x <- abbreviate(colnames(X), abbr.length+6, strict = TRUE, dot = FALSE, method = "both.sides")
  # names_zu <-  abbreviate(colnames(Zu), 9, strict = TRUE, dot = FALSE)
  # names_zv <-  abbreviate(colnames(Zv), 9, strict = TRUE, dot = FALSE)
  # Zu_colnames <- abbreviate(colnames(Zu), abbr.length, strict = TRUE, dot = FALSE, method = "both.sides")
  names_zv <- paste0("lnVARv0i_", abbreviate(colnames(zsv), abbr.length, strict = TRUE, dot = FALSE, method = "both.sides"))


  coef.names.full <- c(
    names_x,
    paste("lnVARv0i_",c("(Intercept)", names_zv[-1]),"", sep = "")
  )
  # cat.print(coef.names.full)
  coef.names.full <- c(
    names_x,  names_zv
  )

  # starting values ---------------------------------------------------------
  
  tymch2 <- lm(Y ~ X - 1)
  olsres  <- resid(tymch2)

  if(ksv==1){
    sv.initial <- log(mean(olsres^2))
  } else {
    sv.initial <- c(log(mean(olsres^2)), rep(0,ksv-1))
  }
  
  theta0 <- c( coef(tymch2), sv.initial)
  names(theta0) <- coef.names.full
  
  if(print.level >= 1){ cat.print(theta0) }
  
  if(print.level >= 1){
    max.name.length <- max(nchar(names(theta0) ))
    est.rez.left <- floor( (max.name.length+42-33) / 2 )
    est.rez.right <- max.name.length+42-33 - est.rez.left
    cat("\n",rep("-", est.rez.left)," Regression with normal errors: ",rep("-", est.rez.right),"\n\n", sep ="")
  }
  # mytrace <- ifelse(print.level < 2, 0, 1)
  if(print.level <= 2){
    trace <- trace1 <- 0
  } else {
    trace1 <- 10
  }
  tymch <- optim(
    par = theta0, fn = .ll.lm.mle, gr = .gr.lm.mle, y = Y, x = X, zsv = zsv, k = k, ksv = ksv,
    method = technique, control = list(fnscale = -1, trace = trace1, maxit = 10000), 
    hessian = TRUE)
  
  tymch$coef <- tymch$par
  tymch$vcov <- tryCatch(solve(-tymch$hessian), tol = .Machine$double.xmin * 10, error = function(e) e )
  
  
  # cat.print(tymch1)
  # table with coefs
  # tymch$par  <- c(tymch$par, sv = unname(sqrt(exp(tymch$par[k+1]))) )
  # tymch$sds  <- suppressWarnings( sqrt(diag(solve(-tymch$hessian))) )
  # tymch$sds  <- c(tymch$sds, tymch$sds[k+1]*exp(tymch$par[k+1]))
  
  tymch$sds  <- suppressWarnings( sqrt(diag(solve(-tymch$hessian))) )

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
  # tymch $ RSS     <- crossprod(eps0)
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
  # tymch$ skewness   <- as.vector(unname(zsk %*%      tymch$par[(k + ksv + 1)            :(k + ksv + ksk)]))
  # tymch$ u     <- as.vector(unname(u))
  # tymch$ eff   <- exp(-tymch$ u)
  
  # cat.print(tymch$coef)
  # cat.print(tymch$par)
  # cat.print(tymch$vcov)
  # cat.print(coef.names.full)
  
  names(tymch$coef) <- names(tymch$par) <- colnames(tymch$vcov) <- rownames(tymch$vcov) <-
    coef.names.full
  
  class(tymch) <- "snnoise"
  return(tymch)
  
}




#


