snreg <- function (
  formula, data, subset,
  start.sk  = NULL,#.5*(2*prod-1),
  ln.var.v  = NULL,
  skew.v    = NULL,
  start.val = NULL,
  technique = c('nr'),#,'bhhh','nm', 'bfgs', 'cg'),
  vcetype   = c('aim'),#, 'opg'), # `approximated information matrix` or `outer product of gradients`
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


