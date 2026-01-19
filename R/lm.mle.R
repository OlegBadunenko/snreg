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


