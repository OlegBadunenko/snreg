
# * T Owen grad -----------------------------------------------------------

TOwen.gr.2 <- function(x,y){
  0.5/pi/(1+y^2.0)*exp(-0.5*x^2.0*(1.0+y^2))
}

TOwen.gr.1 <- function(x,y){
  # -0.5/sqrt(2.0*pi)*(erf(x*y/sqrt(2.0)))*exp(-0.5*x^2.0)
  -0.5/sqrt(2.0*pi)*(2.0*pnorm(x*y)-1.0)*exp(-0.5*x^2.0)
}


# * lm mle ----------------------------------------------------------------

lmd <- function(z) dnorm(z)/pnorm(z)

.ll.lm.mle <- function(theta, y, x, k, zsv, ksv){
  eps0  <- y - x %*% theta[1:k]
  sv    <- sqrt(exp(zsv %*% theta[(k+1):(k+ksv)]))
  # sv2    <- exp(zsv %*% theta[(k+1):(k+ksv)])
  # cat.print(cbind(eps0,sv,al0,epsr,dnorm(epsr/sv, log = TRUE),pnorm(al0*epsr/sv, log.p = TRUE)))
  # ll0   <- -0.5 * ( log(2*pi*sv2) + eps0^2/sv2 )
  ll0   <- -log(sv) + dnorm(eps0/sv, log = TRUE)
  return(sum(ll0))
}

.gr.lm.mle <- function(theta, y, x, k, zsv, ksv){
  eps0  <- y - x %*% theta[1:k]
  sv    <- sqrt(exp(zsv %*% theta[(k+1):(k+ksv)]))
  gr <- cbind(
    sweep(x, 1, eps0/sv^2, FUN = "*"),
    sweep(-0.5*zsv, 1, 1 - eps0^2/sv^2, FUN = "*")
  )
  colnames(gr) <- names(theta)
  colSums( gr )
}

ll.n.gr.b <- function(b, gv, y, x, zv){
  eps0  <- y - x * b
  sv    <- sqrt(exp(zv * gv))
  x*(eps0/sv^2)
}

ll.n.gr.gv <- function(b, gv, y, x, zv){
  eps0  <- y - x * b
  sv    <- sqrt(exp(zv * gv))
  -0.5*zv*(1 - eps0^2/sv^2)
}

# * SN noise --------------------------------------------------------------

.ll.sn0 <- function(theta, y, x, k, ...){
  eps0  <- y - x %*% theta[1:k]
  sv    <- sqrt(exp(theta[k+1]))
  al0   <- theta[k+2]
  epsr  <- eps0 + sv * sqrt(2/pi) * al0 / sqrt(1+al0^2)
  # cat.print(cbind(eps0,sv,al0,epsr,dnorm(epsr/sv, log = TRUE),pnorm(al0*epsr/sv, log.p = TRUE)))
  ll0   <- log(2) - log(sv) + dnorm(epsr/sv, log = TRUE) + pnorm(al0*epsr/sv, log.p = TRUE)
  return(sum(ll0))
}

.ll.sn <- function(theta, y, x, zsv, zsk, k, ksv, ksk, ...){
  eps0  <- y - x %*% theta[1:k]
  # cat.print(eps0)
  sv    <- sqrt(exp(zsv %*% theta[(k+1):(k+ksv)]))
  # cat.print(sv)
  al0   <- zsk %*% theta[(k+ksv+1):(k+ksv+ksk)]
  # cat.print(al0)
  epsr  <- eps0 + sv * sqrt(2/pi) * al0 / sqrt(1+al0^2)
  # cat.print(epsr)
  # cat.print(cbind(eps0,sv,al0,epsr,dnorm(epsr/sv, log = TRUE),pnorm(al0*epsr/sv, log.p = TRUE)))
  ll0   <- log(2) - log(sv) + dnorm(epsr/sv, log = TRUE) + pnorm(al0*epsr/sv, log.p = TRUE)
  # cat.print(ll0)
  return(sum(ll0))
}

.gr.sn <- function(theta, y, x, zsv, zsk, k, ksv, ksk, ...){
  eps0  <- y - x %*% theta[1:k]
  sv    <- sqrt(exp(zsv %*% theta[(k+1):(k+ksv)]))
  al0   <- zsk %*% theta[(k+ksv+1):(k+ksv+ksk)]
  epsr  <- eps0 + sv * sqrt(2/pi) * al0 / sqrt(1+al0^2)
  gr <- cbind(
    sweep(x, 1, epsr/sv^2 - al0/sv*lmd(al0*epsr/sv), FUN = "*"),
    sweep(-0.5*zsv, 1, 1 - epsr/sv^2 * eps0 + lmd(al0*epsr/sv) * al0*eps0/sv, FUN = "*"),
    sweep(zsk, 1, -epsr/ sv* sqrt(2/pi)/ sqrt(1+al0^2)/ (1+al0^2) + lmd(al0*epsr/sv) *( eps0/sv + sqrt(2/pi) * al0 / sqrt(1+al0^2)* (2+al0^2)/(1+al0^2) ), FUN = "*")
  )
  colnames(gr) <- names(theta)
  colSums( gr )
}

.gr.sn.by.i <- function(theta, y, x, zsv, zsk, k, ksv, ksk, ...){
  eps0  <- y - x %*% theta[1:k]
  sv    <- sqrt(exp(zsv %*% theta[(k+1):(k+ksv)]))
  al0   <- zsk %*% theta[(k+ksv+1):(k+ksv+ksk)]
  epsr  <- eps0 + sv * sqrt(2/pi) * al0 / sqrt(1+al0^2)
  gr <- cbind(
    sweep(x, 1, epsr/sv^2 - al0/sv*lmd(al0*epsr/sv), FUN = "*"),
    sweep(-0.5*zsv, 1, 1 - epsr/sv^2 * eps0 + lmd(al0*epsr/sv) * al0*eps0/sv, FUN = "*"),
    sweep(zsk, 1, -epsr/ sv* sqrt(2/pi)/ sqrt(1+al0^2)/ (1+al0^2) + lmd(al0*epsr/sv) *( eps0/sv + sqrt(2/pi) * al0 / sqrt(1+al0^2)* (2+al0^2)/(1+al0^2) ), FUN = "*")
  )
  colnames(gr) <- names(theta)
  return( gr )
}


ll.sn.gr.b <- function(b, gv, ga, y, x, zv, za){
  eps0  <- y - x * b
  sv    <- sqrt(exp(zv * gv))
  al0   <- za * ga
  epsr  <- eps0 + sv * sqrt(2/pi) * al0 / sqrt(1+al0^2)
  x*(epsr/sv^2 - al0/sv*lmd(al0*epsr/sv))
}

ll.sn.gr.gv <- function(b, gv, ga, y, x, zv, za){
  eps0  <- y - x * b
  sv    <- sqrt(exp(zv * gv))
  al0   <- za * ga
  epsr  <- eps0 + sv * sqrt(2/pi) * al0 / sqrt(1+al0^2)
  -0.5*zv*(1 - epsr/sv^2 * eps0 + lmd(al0*epsr/sv) * al0*eps0/sv)
}

ll.sn.gr.ga <- function(b, gv, ga, y, x, zv, za){
  eps0  <- y - x * b
  sv    <- sqrt(exp(zv * gv))
  al0   <- za * ga
  epsr  <- eps0 + sv * sqrt(2/pi) * al0 / sqrt(1+al0^2)
  za*(-epsr/ sv* sqrt(2/pi)/ sqrt(1+al0^2)/ (1+al0^2) + lmd(al0*epsr/sv) *( eps0/sv + sqrt(2/pi) * al0 / sqrt(1+al0^2)* (2+al0^2)/(1+al0^2) )
  )
}

.hess.sn <- function(theta, y, x, zsv, zsk, k, ksv, ksk, ...){
  
  h0 <-  hessian2(funcg = .gr.sn, at = theta, y = y, x = x, zsv=zsv, zsk=zsk, k=k, ksv=ksv, ksk=ksk) 
  
  ifelse(is.na(h0), 0, h0)
  
}

.hess.gr.sn <- function(theta, y, x, zsv, zsk, k, ksv, ksk, ...){
  
  g0 <- .gr.sn(theta, y = y, x = x, zsv=zsv, zsk=zsk, k=k, ksv=ksv, ksk=ksk)
  
  h0 <-  hessian2(funcg = .gr.sn, at = theta, y = y, x = x, zsv=zsv, zsk=zsk, k=k, ksv=ksv, ksk=ksk) 
  
  list( grad = ifelse(is.na(g0), 0, g0),  hessian1 = ifelse(is.na(h0), 0, h0) )
  
}


# SN-Exp ------------------------------------------------------------------

# log likelihood ----------------------------------------------------------

.ll.sn.exp.R <- function(theta, myprod, y, x, zsv, zsk, zsu, n.ids, k, ksv, ksk, ksu, threads = 1, ...){
  # k:    betas
  # ksv:  noise, scale/variance
  # ksk:  noise, skew
  # ksu:  inefficiency: scale/variance
  # cat.print(myprod)
  # cat.print(1)
  
  eps0  <- y - x %*%    theta[1                        :k]
  # cat.print(2)
  sv    <- sqrt(exp( zsv %*% theta[(k + 1)                  :(k + ksv)]))
  # cat.print(sv2)
  al0   <- zsk %*%      theta[(k + ksv + 1)            :(k + ksv + ksk)]
  # cat.print(4)
  lam   <- 1 / sqrt(exp( zsu %*% theta[(k + ksv + ksk  + 1)     :(k + ksv + ksk + ksu)] ))
  # cat.print(su2)
  # cat.print(6)
  
  epsr  <- eps0 + sv * sqrt(2/pi) * al0 / sqrt(1 + al0^2)
  # print(summary(epsr))
  
  u1    <- (myprod*epsr + lam*sv^2)/sv
  b1    <- al0 * myprod
  a1    <- -b1 * lam * sv
  a2    <- a1 / sqrt(1+b1^2)
  # a2    <- ifelse(a1 == 0, 1, a2_)
  # mupr  <- -mu1r/sstar
  
  # ll    <- numeric(n)
  # ll0   <- 0.5*( pnorm(a2) - pnorm(u1) + 1*(a2/u1<0) ) #terms 1, 2, and 5
  
  # cat.print(11)
  # as    <- cbind(a1 /mupr + b1, mupr * a1 / a2^2 + b1)
  # # as    <- ifelse(is.na(as), 0, as)
  # if(sum(is.na(as)) > 0){
  #   tymch1 <- data.frame(eps0=unname(eps0),sv2,al0,su2,mu0,epsr=unname(epsr),a1,b1,mupr=unname(mupr),a2=unname(a2),ll0 = unname(ll0))
  #   cat.print(head(tymch1,17))
  #   cat.print(theta)
  #   warning("Something is wrong, NA")
  #   return(tymch1)
  # }
  # cat("Count infinite ",sum(is.infinite(as)),"; Count NA ",sum(is.na(as)),"\n")
  # cat.print(dim(as))
  
  # ll1   <- ll0
  # term1 <- term2 <- term3 <- term4 <- numeric(n.ids)
  # for (i in 1:n.ids) {
  #   term1[i] <- -sn::T.Owen( u1[i], a2[i] / u1[i] )
  #   term2[i] <- -sn::T.Owen( a2[i], u1[i] / a2[i] )
  #   term3[i] <- sn::T.Owen( u1[i], b1[i] + a1[i] / u1[i] )
  #   term4[i] <- sn::T.Owen( a2[i], b1[i] + u1[i] * (1 + b1[i]^2) / a1[i] )
  # }
  
  
  # toInclude <- is.finite(u1) & is.finite(a2) & is.finite(a2)
  # n <- length(y)
  # # u1_na <- is.na(u1) | is.infinite(u1)
  # # a2_na <- is.na(a2) | is.infinite(a2)
  # if(sum(toInclude) < n){
  #   # toInclude <- !(u1_na + a2_na)
  #   term1 <- term2 <- term3 <- term4 <- rep(.Machine$double.xmin, n)
  #   term1[toInclude] <- -TOwen( u1[toInclude], a2[toInclude] / u1[toInclude], threads )    
  #   term2[toInclude] <- -TOwen( a2[toInclude], u1[toInclude] / a2[toInclude], threads )
  #   term3[toInclude] <- TOwen( u1[toInclude], b1[toInclude] + a1[toInclude] / u1[toInclude], threads )
  #   term4[toInclude] <- TOwen( a2[toInclude], b1[toInclude] + u1[toInclude] * (1 + b1[toInclude]^2) / a1[toInclude], threads )
  # } else {
  #   term1 <- -TOwen( u1, a2 / u1, threads )
  #   term2 <- -TOwen( a2, u1 / a2, threads )
  #   term3 <- TOwen( u1, b1 + a1 / u1, threads )
  #   term4 <- TOwen( a2, b1 + u1 * (1 + b1^2) / a1, threads )
  #   # term5 <- pnorm(a2) * pnorm(-u1)
  # }
  
  term1 <- -TOwen( u1, a2 / u1, threads )
  term2 <- -TOwen( a2, u1 / a2, threads )
  term3 <- TOwen( u1, b1 + a1 / u1, threads )
  term4 <- TOwen( a2, b1 + u1 * (1 + b1^2) / a1, threads )
  term5 <- pnorm(a2) * pnorm(-u1)
  t1_5  <- term1 + term2 + term3 + term4 + term5
  # cat.print(any(is.na(t1_5)))
  # tymch <- cbind(term1 , term2 , term3 , term4, term5, t1_5)
  # cat.print(tymch[1:20,])
  
  # t1_5  <- ifelse(t1_5 <= 0 | !is.finite(t1_5), .Machine$double.xmin, t1_5)
  t1_5  <- ifelse(t1_5 <= 0, .Machine$double.xmin, t1_5)
  # tymch <- cbind(term1 , term2 , term3 , term4, term5, t1_5)
  # cat.print(tymch[1:20,])
  # cat.print(t(t1_5))
  
  # if(any(t1_5 < 0)) cat.print(as.vector(t1_5))
  # t1_5 <- ifelse(t1_5 < 0, .Machine$double.eps, t1_5 )
  # print[(head(cbind(ll0, term3, term4, ll0 + term3 + term4)))
  # tymch1 <- data.frame(eps0=unname(eps0),epsr=unname(epsr),sv,lam,al0,a1,a2=unname(a2),ll0 = unname(ll0),
  #                      ll1 = unname(ll0 + term3 + term4),
  #                      terms0 = unname(log(2*lam) + myprod*epsr*lam + lam^2*sv^2/2)
  # )
  # cat.print(head(tymch1,3))
  # print(head(cbind(ll0, term3, term4, ll0 + term3 + term4)))
  # ll00 <- ll0
  # ll0 <- ll0 + term3 + term4
  # ll1 <- ifelse(ll0 < 0 | abs(ll0) < .Machine$double.eps, .Machine$double.eps, ll0)
  # cat.print(t(ll))
  # log of product begin
  # ll2    <- ifelse(abs(ll1) < 1e-16 | ll1 < 0, .Machine$double.eps, ll1)
  # cat.print(t(ll))
  # tymch1 <- 2 / pnorm(mu0/sqrt(su2)) * sstar / sqrt(s2) *
  #               dnorm( (mu0 + myprod * epsr) / sqrt(s2) )
  # tymch2 <- ifelse(abs(tymch1) < 1e-16 | tymch1 < 0, .Machine$double.eps, tymch1)
  # ll3    <- log( ll2 * tymch2 )
  # ll     <- sum(ll3)
  # # log of product end
  # if(is.infinite(ll)){
  #   tymch1 <- data.frame(eps0=unname(eps0),sv2,al0,su2,mu0,
  #                        epsr=unname(epsr),a1,b1,mupr=unname(mupr),a2,
  #                        tymch1=unname(tymch1),ll0=unname(ll0),
  #                        ll1=unname(ll1),ll2=unname(ll2),ll3=unname(ll3))
  #   cat.print(head(tymch1),30)
  #   cat.print(theta)
  #   warning("Something is wrong, Inf")
  #   return(tymch1)
  # }
  # sum of logs begin
  # print(table(ll0 < 0))
  # tymch3 <- cbind(log( t1_5 ) , log(2*lam) , myprod*epsr*lam + lam^2*sv^2/2)
  # cat.print(tymch3[1:20,])
  # tymch2 <- is.finite(tymch3)
  # toInclude <- rowSums(tymch2) < 3
  # cat.print(tymch3[toInclude,])
  # cat.print(tymch[toInclude,])
  tymch1 <- log( t1_5 ) + log(2*lam) + myprod*epsr*lam + lam^2*sv^2/2
  # cat.print(sum(tymch1))
  return( sum(tymch1, na.rm = TRUE) )
  # ll4    <- log( (term1 + term2 + term3 + term4 + term5) * 2*lam * exp(myprod*epsr*lam + lam^2*sv^2/2))
  # tymch2 <- data.frame(eps0=unname(eps0),epsr=unname(epsr),
  #                      sv,lam,
  #                      al0,a1,a2=unname(a2),
  #                      bau1 = unname(b1 + a1 / u1),
  #                      bub2a = unname(b1 + u1 * (1 + b1^2) / a1),
  #                      u1 = unname(u1),
  #                      t125 =  unname(ll00),
  #                      t3 = unname(term3),
  #                      t4 = unname(term4),
  #                      t1_5before = unname(ll0),
  #                      t1_5before1 = unname(term1 + term2 + term3 + term4 + term5),
  #                      t1_5after  = unname(ll1),
  #                      llpar2 = unname(log(2*lam) + myprod*epsr*lam + lam^2*sv^2/2),
  #                      ll3 = unname(ll3), 
  #                      ll4 = unname(ll4)
  # )
  # cat.print(head(tymch2,20))
  # print(table(is.infinite(ll3)))
  # ll     <- sum(ll3)
  # cat.print(as.vector(ll3))
  # sum of logs end
  # return( tymch2)
  # cat.print(ll3)
  # ll3 <- ifelse(is.infinite(ll3) & ll3 < 0, NaN, ll3 )
  # ll <- sum(tymch1, na.rm = TRUE)
  # if(ll4 > 5000){
  #   t1_5.pos <- term1 + term2 + term3 + term4 + term5 > 0
  #   t1_5.min.pos <- min((term1 + term2 + term3 + term4 + term5)[t1_5.pos])
  #   tymch2 <- cbind(t1_5.pos, term1 + term2 + term3 + term4 + term5, t1_5, log( t1_5 ), log(2*lam) + myprod*epsr*lam + lam^2*sv^2/2, ll3)
  #   colnames(tymch2) <- c("t1_5.pos", "t1_5_0", "t1_5_1", "lnINT", "ll0", "ll")
  #   cat.print(tymch2[1:51,])
  #   cat.print( min((term1 + term2 + term3 + term4 + term5)[t1_5.pos]) )
  #   cat.print(ll3[t1_5.pos])
  #   cat.print(ll3[!t1_5.pos])
  #   Sys.sleep(5)
  # } 
  # cat.print(ll3)
  # cat.print(sum(ll3))
  # return( ll )
}

.ll.sn.exp.R.by.i <- function(theta, myprod, y, x, zsv, zsk, zsu, n.ids, k, ksv, ksk, ksu, threads = 1, ...){
  # k:    betas
  # ksv:  noise, scale/variance
  # ksk:  noise, skew
  # ksu:  inefficiency: scale/variance
  # cat.print(myprod)
  # cat.print(1)
  
  eps0  <- y - x %*%    theta[1                        :k]
  # cat.print(2)
  sv    <- sqrt(exp( zsv %*% theta[(k + 1)                  :(k + ksv)]))
  # cat.print(sv2)
  al0   <- zsk %*%      theta[(k + ksv + 1)            :(k + ksv + ksk)]
  # cat.print(4)
  lam   <- 1 / sqrt(exp( zsu %*% theta[(k + ksv + ksk  + 1)     :(k + ksv + ksk + ksu)] ))
  # cat.print(su2)
  # cat.print(6)
  
  epsr  <- eps0 + sv * sqrt(2/pi) * al0 / sqrt(1 + al0^2)
  # print(summary(epsr))
  
  u1    <- (myprod*epsr + lam*sv^2)/sv
  b1    <- al0 * myprod
  a1    <- -b1 * lam * sv
  a2    <- a1 / sqrt(1+b1^2)
  
  # tymch <- cbind(u1,a1,a2,b1)
  # cat.print(tymch[1:10,])
  # a2    <- ifelse(a1 == 0, 1, a2_)
  # mupr  <- -mu1r/sstar
  
  # ll    <- numeric(n)
  # ll0   <- 0.5*( pnorm(a2) - pnorm(u1) + 1*(a2/u1<0) ) #terms 1, 2, and 5
  
  # cat.print(11)
  # as    <- cbind(a1 /mupr + b1, mupr * a1 / a2^2 + b1)
  # # as    <- ifelse(is.na(as), 0, as)
  # if(sum(is.na(as)) > 0){
  #   tymch1 <- data.frame(eps0=unname(eps0),sv2,al0,su2,mu0,epsr=unname(epsr),a1,b1,mupr=unname(mupr),a2=unname(a2),ll0 = unname(ll0))
  #   cat.print(head(tymch1,17))
  #   cat.print(theta)
  #   warning("Something is wrong, NA")
  #   return(tymch1)
  # }
  # cat("Count infinite ",sum(is.infinite(as)),"; Count NA ",sum(is.na(as)),"\n")
  # cat.print(dim(as))
  
  # ll1   <- ll0
  # term1 <- term2 <- term3 <- term4 <- numeric(n.ids)
  # for (i in 1:n.ids) {
  #   term1[i] <- -sn::T.Owen( u1[i], a2[i] / u1[i] )
  #   term2[i] <- -sn::T.Owen( a2[i], u1[i] / a2[i] )
  #   term3[i] <- sn::T.Owen( u1[i], b1[i] + a1[i] / u1[i] )
  #   term4[i] <- sn::T.Owen( a2[i], b1[i] + u1[i] * (1 + b1[i]^2) / a1[i] )
  # }
  toInclude <- is.finite(u1) & is.finite(a2) & is.finite(a2)
  n <- length(y)
  # u1_na <- is.na(u1) | is.infinite(u1)
  # a2_na <- is.na(a2) | is.infinite(a2)
  if(sum(toInclude) < n){
    # toInclude <- !(u1_na + a2_na)
    term1 <- term2 <- term3 <- term4 <- rep(.Machine$double.xmin, n)
    term1[toInclude] <- -TOwen( u1[toInclude], a2[toInclude] / u1[toInclude], threads )    
    term2[toInclude] <- -TOwen( a2[toInclude], u1[toInclude] / a2[toInclude], threads )
    term3[toInclude] <- TOwen( u1[toInclude], b1[toInclude] + a1[toInclude] / u1[toInclude], threads )
    term4[toInclude] <- TOwen( a2[toInclude], b1[toInclude] + u1[toInclude] * (1 + b1[toInclude]^2) / a1[toInclude], threads )
  } else {
    term1 <- -TOwen( u1, a2 / u1, threads )
    term2 <- -TOwen( a2, u1 / a2, threads )
    term3 <- TOwen( u1, b1 + a1 / u1, threads )
    term4 <- TOwen( a2, b1 + u1 * (1 + b1^2) / a1, threads )
    # term5 <- pnorm(a2) * pnorm(-u1)
  }
  term5 <- pnorm(a2) * pnorm(-u1)
  t1_5  <- term1 + term2 + term3 + term4 + term5
  
  t1_5  <- ifelse(t1_5 < 0, .Machine$double.xmin, t1_5)
  
  # if(any(t1_5 < 0)) cat.print(as.vector(t1_5))
  # t1_5 <- ifelse(t1_5 < 0, .Machine$double.eps, t1_5 )
  # print[(head(cbind(ll0, term3, term4, ll0 + term3 + term4)))
  # tymch1 <- data.frame(eps0=unname(eps0),epsr=unname(epsr),sv,lam,al0,a1,a2=unname(a2),ll0 = unname(ll0),
  #                      ll1 = unname(ll0 + term3 + term4),
  #                      terms0 = unname(log(2*lam) + myprod*epsr*lam + lam^2*sv^2/2)
  # )
  # cat.print(head(tymch1,3))
  # print(head(cbind(ll0, term3, term4, ll0 + term3 + term4)))
  # ll00 <- ll0
  # ll0 <- ll0 + term3 + term4
  # ll1 <- ifelse(ll0 < 0 | abs(ll0) < .Machine$double.eps, .Machine$double.eps, ll0)
  # cat.print(t(ll))
  # log of product begin
  # ll2    <- ifelse(abs(ll1) < 1e-16 | ll1 < 0, .Machine$double.eps, ll1)
  # cat.print(t(ll))
  # tymch1 <- 2 / pnorm(mu0/sqrt(su2)) * sstar / sqrt(s2) *
  #               dnorm( (mu0 + myprod * epsr) / sqrt(s2) )
  # tymch2 <- ifelse(abs(tymch1) < 1e-16 | tymch1 < 0, .Machine$double.eps, tymch1)
  # ll3    <- log( ll2 * tymch2 )
  # ll     <- sum(ll3)
  # # log of product end
  # if(is.infinite(ll)){
  #   tymch1 <- data.frame(eps0=unname(eps0),sv2,al0,su2,mu0,
  #                        epsr=unname(epsr),a1,b1,mupr=unname(mupr),a2,
  #                        tymch1=unname(tymch1),ll0=unname(ll0),
  #                        ll1=unname(ll1),ll2=unname(ll2),ll3=unname(ll3))
  #   cat.print(head(tymch1),30)
  #   cat.print(theta)
  #   warning("Something is wrong, Inf")
  #   return(tymch1)
  # }
  # sum of logs begin
  # print(table(ll0 < 0))
  tymch1 <- log( t1_5 ) + log(2*lam) + myprod*epsr*lam + lam^2*sv^2/2
  return( tymch1 )
  # ll4    <- log( (term1 + term2 + term3 + term4 + term5) * 2*lam * exp(myprod*epsr*lam + lam^2*sv^2/2))
  # tymch2 <- data.frame(eps0=unname(eps0),epsr=unname(epsr),
  #                      sv,lam,
  #                      al0,a1,a2=unname(a2),
  #                      bau1 = unname(b1 + a1 / u1),
  #                      bub2a = unname(b1 + u1 * (1 + b1^2) / a1),
  #                      u1 = unname(u1),
  #                      t125 =  unname(ll00),
  #                      t3 = unname(term3),
  #                      t4 = unname(term4),
  #                      t1_5before = unname(ll0),
  #                      t1_5before1 = unname(term1 + term2 + term3 + term4 + term5),
  #                      t1_5after  = unname(ll1),
  #                      llpar2 = unname(log(2*lam) + myprod*epsr*lam + lam^2*sv^2/2),
  #                      ll3 = unname(ll3), 
  #                      ll4 = unname(ll4)
  # )
  # cat.print(head(tymch2,20))
  # print(table(is.infinite(ll3)))
  # ll     <- sum(ll3)
  # cat.print(as.vector(ll3))
  # sum of logs end
  # return( tymch2)
  # cat.print(ll3)
  # ll3 <- ifelse(is.infinite(ll3) & ll3 < 0, NaN, ll3 )
  # ll <- sum(tymch1, na.rm = TRUE)
  # if(ll4 > 5000){
  #   t1_5.pos <- term1 + term2 + term3 + term4 + term5 > 0
  #   t1_5.min.pos <- min((term1 + term2 + term3 + term4 + term5)[t1_5.pos])
  #   tymch2 <- cbind(t1_5.pos, term1 + term2 + term3 + term4 + term5, t1_5, log( t1_5 ), log(2*lam) + myprod*epsr*lam + lam^2*sv^2/2, ll3)
  #   colnames(tymch2) <- c("t1_5.pos", "t1_5_0", "t1_5_1", "lnINT", "ll0", "ll")
  #   cat.print(tymch2[1:51,])
  #   cat.print( min((term1 + term2 + term3 + term4 + term5)[t1_5.pos]) )
  #   cat.print(ll3[t1_5.pos])
  #   cat.print(ll3[!t1_5.pos])
  #   Sys.sleep(5)
  # } 
  # cat.print(ll3)
  # cat.print(sum(ll3))
  # return( ll )
}

.ll.sn.exp <- function(theta, myprod, y, x, zsv, zsk, zsu, n.ids, k, ksv, ksk, ksu, threads = 1, ...){
  # cat.print(theta)
  tymch <- .C("ll_sn_exp", as.integer(threads),
              as.double(myprod), as.double(y), as.double(x),
              as.double(zsv), as.double(zsk), as.double(zsu),
              as.integer(n.ids), as.integer(k), 
              as.integer(ksv), as.integer(ksk), as.integer(ksu),
              as.double(theta), 
              lnls = double(n.ids), lnl = as.double(1)) 
  # cat.print(sum(tymch$lnls))
  return(tymch$lnl)
}

.u.sn.exp <- function(theta, myprod, y, x, zsv, zsk, zsu, n.ids, k, ksv, ksk, ksu, threads = 1, ...){
  # cat.print(theta)
  tymch <- .C("u_sn_exp", as.integer(threads),
              as.double(myprod), as.double(y), as.double(x),
              as.double(zsv), as.double(zsk), as.double(zsu),
              as.integer(n.ids), as.integer(k), 
              as.integer(ksv), as.integer(ksk), as.integer(ksu),
              as.double(theta), 
              u = double(n.ids), 
              eps = double(n.ids), 
              sv = double(n.ids), 
              al = double(n.ids), 
              lam = double(n.ids))
  # cat.print(sum(tymch$lnls))
  return(tymch)
}


.ll.sn.exp.by.i <- function(theta, myprod, y, x, zsv, zsk, zsu, n.ids, k, ksv, ksk, ksu, threads = 1, ...){
  tymch <- .C("ll_sn_exp", as.integer(threads),
              as.double(myprod), 
              as.double(y), as.double(x),
              as.double(zsv), as.double(zsk), as.double(zsu),
              as.integer(n.ids), as.integer(k), 
              as.integer(ksv), as.integer(ksk), as.integer(ksu),
              as.double(theta), 
              lnls = double(n.ids), lnl = as.double(1)) 
  # cat.print(tymch$lnl)
  return(tymch$lnls)
}

# gradient ----------------------------------------------------------------

.gr.sn.exp.appr <- function(theta, myprod, y, x, zsv, zsk, zsu, n.ids, k, ksv, ksk, ksu, threads = 1, ...) {
  
  g0 <- gradient1(func = .ll.sn.exp, at = theta, myprod=myprod, y=y, x=x, zsv=zsv, zsk=zsk, zsu=zsu, n.ids=n.ids, k=k, ksv=ksv, ksk=ksk, ksu=ksu, threads = threads)
  names(g0) <- names(theta)
  return(ifelse(is.na(g0) | is.infinite(g0), 0, g0))
  
}

.gr.sn.exp <- function(theta, myprod, y, x, zsv, zsk, zsu, n.ids, k, ksv, ksk, ksu, threads = 1, ...){
  # cat.print(theta)
  tymch <- .C("sn_exp_ll_grad", as.integer(threads),
              as.double(myprod), as.double(y), as.double(x),
              as.double(zsv), as.double(zsk), as.double(zsu),
              as.integer(n.ids), as.integer(k), 
              as.integer(ksv), as.integer(ksk), as.integer(ksu),
              as.double(theta), 
              lnls = double(n.ids), lnl = as.double(1), 
              grad = double(n.ids*(k+ksv+ksk+ksu)), grad2 = double(k+ksv+ksk+ksu), t1_5 = double(n.ids)) 
  # cat.print(sum(tymch$lnls))
  return(tymch$grad2)
}

.gr.sn.exp.R <- function(theta, myprod, y, x, zsv, zsk, zsu, n.ids, k, ksv, ksk, ksu, threads = 1, ...) {
  
  eps0  <- y - x %*% theta[1 :k]
  sv    <- sqrt(exp( zsv %*% theta[(k + 1) :(k + ksv)]))
  alp   <- zsk %*% theta[(k + ksv + 1) :(k + ksv + ksk)]
  lam   <- 1 / sqrt(exp( zsu %*% theta[(k + ksv + ksk  + 1) :(k + ksv + ksk + ksu)] ))
  su    <- 1/ lam
  
  epsr  <- eps0 + sv * sqrt(2/pi) * alp / sqrt(1 + alp^2)
  u1    <- (myprod*epsr + lam*sv^2)/sv
  b     <- alp * myprod
  a     <- -b * lam * sv
  a2    <- a / sqrt(1 + b^2)
  
  term1 <- -TOwen( u1, a2 / u1, threads = 1 )
  term2 <- -TOwen( a2, u1 / a2, threads = 1 )
  term3 <- TOwen( u1, b + a / u1, threads = 1 )
  term4 <- TOwen( a2, b + u1 * (1 + b^2) / a, threads = 1 )
  term5 <- pnorm(a2) * pnorm(-u1)
  t1_5  <- term1 + term2 + term3 + term4 + term5
  
  eb    <- -x
  SVU   <- sv * lam
  SVUv  <- sweep(zsv, 1, SVU * 0.5, FUN = "*") #  SVU * 0.5 * zsv
  SVUu  <- sweep(zsu, 1, -SVU * 0.5, FUN = "*") # -SVU * 0.5 * zsu
  A     <- alp / sqrt(1+alp^2)
  As    <- sweep(zsk, 1, (1+alp^2)^(-3/2), FUN = "*") #zsk * (1-A^2) / sqrt(1+alp^2)
  
  d_u1_b  <- u1b <- sweep(eb, 1, myprod / sv, FUN = "*") # myprod/sv*eb
  d_u1_v  <- u1v <- sweep(zsv, 1, myprod * eps0 * -0.5 / sv, FUN = "*") + SVUv # myprod * eps0 * -0.5 * sv * zsv
  d_u1_s  <- u1s <- myprod * sqrt(2/pi) * As
  d_u1_u  <- u1u <- SVUu
  
  d_a2_b  <- matrix(0, nrow = nrow(x), ncol = ncol(x))
  d_a2_v  <- a2v <- sweep(SVUv, 1, -myprod * A, FUN = "*")  
  d_a2_s  <- a2s <- sweep(As,   1, -myprod * SVU, FUN = "*")  
  d_a2_u  <- a2u <- sweep(SVUu, 1, -myprod * A, FUN = "*")  
  
  d_a2u1_b <-                                  - sweep(u1b, 1, a2 / u1^2, FUN = "*")
  d_a2u1_v <- sweep(a2v, 1, 1 / u1, FUN = "*") - sweep(u1v, 1, a2 / u1^2, FUN = "*")
  d_a2u1_s <- sweep(a2s, 1, 1 / u1, FUN = "*") - sweep(u1s, 1, a2 / u1^2, FUN = "*")
  d_a2u1_u <- sweep(a2u, 1, 1 / u1, FUN = "*") - sweep(u1u, 1, a2 / u1^2, FUN = "*")
  
  d_u1a2_b <- sweep(u1b, 1, 1 / a2, FUN = "*")
  d_u1a2_v <- sweep(u1v, 1, 1 / a2, FUN = "*") - sweep(a2v, 1, u1 / a2^2, FUN = "*")
  d_u1a2_s <- sweep(u1s, 1, 1 / a2, FUN = "*") - sweep(a2s, 1, u1 / a2^2, FUN = "*")
  d_u1a2_u <- sweep(u1u, 1, 1 / a2, FUN = "*") - sweep(a2u, 1, u1 / a2^2, FUN = "*")
  
  av       <- sweep(SVUv, 1, -myprod * alp, FUN = "*")
  alps     <- zsk
  as       <- sweep(alps, 1, -myprod * SVU, FUN = "*")
  au       <- sweep(SVUu, 1, -myprod * alp, FUN = "*")
  
  d_bau1_b <-                                 - sweep(u1b, 1, a / u1^2, FUN = "*")
  d_bau1_v <- sweep(av, 1, 1 / u1, FUN = "*") - sweep(u1v, 1, a / u1^2, FUN = "*")
  d_bau1_s <- sweep(as, 1, 1 / u1, FUN = "*") - sweep(u1s, 1, a / u1^2, FUN = "*") + myprod * alps
  d_bau1_u <- sweep(au, 1, 1 / u1, FUN = "*") - sweep(u1u, 1, a / u1^2, FUN = "*")
  
  d_bu1a_b <- sweep(u1b, 1, (1+b^2) / a, FUN = "*")
  d_bu1a_v <- sweep(u1v, 1, (1+b^2) / a, FUN = "*") - sweep(av, 1, u1*(1+b^2) / a^2, FUN = "*")
  d_bu1a_u <- sweep(u1u, 1, (1+b^2) / a, FUN = "*") - sweep(au, 1, u1*(1+b^2) / a^2, FUN = "*")
  # d_bu1a_s <- sweep(u1s, 1, (1+b^2) / a, FUN = "*") + sweep(alps, 1, u1 * ( a*2*alp + myprod * SVU *(1+b^2) ) / a^2, FUN = "*") + myprod * alps
  # d_bu1a_s <- sweep(u1s, 1, (1+b^2) / a, FUN = "*") + sweep(alps, 1, myprod*2*b*u1/a, FUN = "*") - sweep(as, 1, u1*(1+b^2)/a^2, FUN = "*") + myprod * alps
  # d_bu1a_s <- sweep(u1s, 1, (1+b^2) / a, FUN = "*") - sweep(as, 1, u1*((1+b^2)/a^2 - 2/SVU^2), FUN = "*") + myprod * alps
  # d_bu1a_s <- sweep(u1s, 1, (1+b^2) / a, FUN = "*") - sweep(as, 1, u1*((1+b^2)/a^2 + 2*b/a/SVU), FUN = "*") + myprod * alps
  d_bu1a_s <- sweep(u1s, 1, (1+b^2) / a, FUN = "*") - sweep(as, 1, u1/SVU^2*(1/b^2-1), FUN = "*") + myprod * alps
  
  # d_ln2su_u <- -0.5*zsu
  
  d_rest_b <- sweep(eb,   1, myprod/su, FUN = "*")
  d_rest_v <- sweep(SVUv, 1, myprod*sqrt(2/pi)*A + SVU, FUN = "*")
  d_rest_s <- sweep(As,   1, myprod*sqrt(2/pi)*SVU, FUN = "*")
  d_rest_u <- sweep(SVUu, 1, myprod*sqrt(2/pi)*A + SVU, FUN = "*") + sweep(zsu, 1, myprod * eps0 * -0.5 / su -0.5, FUN = "*") #- 0.5*zsu
  
  T_gr_1_u1_a2u1 <- TOwen.gr.1(u1,a2/u1); T_gr_2_u1_a2u1 <- TOwen.gr.2(u1,a2/u1)
  
  T_gr_1_a2_u1a2 <- TOwen.gr.1(a2,u1/a2); T_gr_2_a2_u1a2 <- TOwen.gr.2(a2,u1/a2)
  
  T_gr_1_u1_bau1 <- TOwen.gr.1(u1,b+a/u1); T_gr_2_u1_bau1 <- TOwen.gr.2(u1,b+a/u1)
  
  T_gr_1_a2_bu1a <- TOwen.gr.1(a2,b+u1*(1+b^2)/a); T_gr_2_a2_bu1a <- TOwen.gr.2(a2,b+u1*(1+b^2)/a)
  
  d_term1_b <- sweep(d_u1_b, 1, T_gr_1_u1_a2u1, FUN = "*") + sweep(d_a2u1_b, 1, T_gr_2_u1_a2u1, FUN = "*") 
  d_term2_b <- sweep(d_a2_b, 1, T_gr_1_a2_u1a2, FUN = "*") + sweep(d_u1a2_b, 1, T_gr_2_a2_u1a2, FUN = "*") 
  d_term3_b <- sweep(d_u1_b, 1, T_gr_1_u1_bau1, FUN = "*") + sweep(d_bau1_b, 1, T_gr_2_u1_bau1, FUN = "*") 
  d_term4_b <- sweep(d_a2_b, 1, T_gr_1_a2_bu1a, FUN = "*") + sweep(d_bu1a_b, 1, T_gr_2_a2_bu1a, FUN = "*")
  
  d_term1_v <- sweep(d_u1_v, 1, T_gr_1_u1_a2u1, FUN = "*") + sweep(d_a2u1_v, 1, T_gr_2_u1_a2u1, FUN = "*") 
  d_term2_v <- sweep(d_a2_v, 1, T_gr_1_a2_u1a2, FUN = "*") + sweep(d_u1a2_v, 1, T_gr_2_a2_u1a2, FUN = "*") 
  d_term3_v <- sweep(d_u1_v, 1, T_gr_1_u1_bau1, FUN = "*") + sweep(d_bau1_v, 1, T_gr_2_u1_bau1, FUN = "*") 
  d_term4_v <- sweep(d_a2_v, 1, T_gr_1_a2_bu1a, FUN = "*") + sweep(d_bu1a_v, 1, T_gr_2_a2_bu1a, FUN = "*")
  
  d_term1_s <- sweep(d_u1_s, 1, T_gr_1_u1_a2u1, FUN = "*") + sweep(d_a2u1_s, 1, T_gr_2_u1_a2u1, FUN = "*") 
  d_term2_s <- sweep(d_a2_s, 1, T_gr_1_a2_u1a2, FUN = "*") + sweep(d_u1a2_s, 1, T_gr_2_a2_u1a2, FUN = "*") 
  d_term3_s <- sweep(d_u1_s, 1, T_gr_1_u1_bau1, FUN = "*") + sweep(d_bau1_s, 1, T_gr_2_u1_bau1, FUN = "*") 
  d_term4_s <- sweep(d_a2_s, 1, T_gr_1_a2_bu1a, FUN = "*") + sweep(d_bu1a_s, 1, T_gr_2_a2_bu1a, FUN = "*")
  
  d_term1_u <- sweep(d_u1_u, 1, T_gr_1_u1_a2u1, FUN = "*") + sweep(d_a2u1_u, 1, T_gr_2_u1_a2u1, FUN = "*") 
  d_term2_u <- sweep(d_a2_u, 1, T_gr_1_a2_u1a2, FUN = "*") + sweep(d_u1a2_u, 1, T_gr_2_a2_u1a2, FUN = "*") 
  d_term3_u <- sweep(d_u1_u, 1, T_gr_1_u1_bau1, FUN = "*") + sweep(d_bau1_u, 1, T_gr_2_u1_bau1, FUN = "*") 
  d_term4_u <- sweep(d_a2_u, 1, T_gr_1_a2_bu1a, FUN = "*") + sweep(d_bu1a_u, 1, T_gr_2_a2_bu1a, FUN = "*")
  
  Phis_gr_1 <- pnorm(a2) * dnorm(-u1)
  Phis_gr_2 <- pnorm(-u1) * dnorm(a2)
  
  d_term5_b <- sweep(-d_u1_b, 1, Phis_gr_1, FUN = "*") + sweep(d_a2_b, 1, Phis_gr_2, FUN = "*")
  d_term5_v <- sweep(-d_u1_v, 1, Phis_gr_1, FUN = "*") + sweep(d_a2_v, 1, Phis_gr_2, FUN = "*")
  d_term5_s <- sweep(-d_u1_s, 1, Phis_gr_1, FUN = "*") + sweep(d_a2_s, 1, Phis_gr_2, FUN = "*")
  d_term5_u <- sweep(-d_u1_u, 1, Phis_gr_1, FUN = "*") + sweep(d_a2_u, 1, Phis_gr_2, FUN = "*")
  
  g_b <- d_rest_b + sweep(-d_term1_b-d_term2_b+d_term3_b+d_term4_b+d_term5_b, 1, t1_5, FUN = "/")
  g_v <- d_rest_v + sweep(-d_term1_v-d_term2_v+d_term3_v+d_term4_v+d_term5_v, 1, t1_5, FUN = "/")
  g_s <- d_rest_s + sweep(-d_term1_s-d_term2_s+d_term3_s+d_term4_s+d_term5_s, 1, t1_5, FUN = "/")
  g_u <- d_rest_u + sweep(-d_term1_u-d_term2_u+d_term3_u+d_term4_u+d_term5_u, 1, t1_5, FUN = "/")
  
  g0  <- colSums(cbind(g_b, g_v, g_s, g_u), na.rm = TRUE)
  
  return(unname(g0))
  
}

.gr.sn.exp.by.i.appr <- function(theta, myprod, y, x, zsv, zsk, zsu, n.ids, k, ksv, ksk, ksu, threads = 1, ...) {
  # cat("this is .gr.sn.exp.by.i  \n\n")
  # cat.print(theta)
  # 
  g0 <- gradient1.by.i(func = .ll.sn.exp.by.i, at = theta, myprod=myprod, y=y, x=x, zsv=zsv, zsk=zsk, zsu=zsu, n.ids=n.ids, k=k, ksv=ksv, ksk=ksk, ksu=ksu, threads = threads)
  # cat("this is .gr.sn.exp.by.i  \n\n")
  # cat.print(g0)
  colnames(g0) <- names(theta)
  return(ifelse(is.na(g0) | is.infinite(g0), 0, g0))
  
}

.gr.sn.exp.by.i <- function(theta, myprod, y, x, zsv, zsk, zsu, n.ids, k, ksv, ksk, ksu, threads = 1, ...){
  # cat.print(theta)
  tymch <- .C("sn_exp_ll_grad", as.integer(threads),
              as.double(myprod), as.double(y), as.double(x),
              as.double(zsv), as.double(zsk), as.double(zsu),
              as.integer(n.ids), as.integer(k), 
              as.integer(ksv), as.integer(ksk), as.integer(ksu),
              as.double(theta), 
              lnls = double(n.ids), lnl = as.double(1), 
              grad = double(n.ids*(k+ksv+ksk+ksu)), grad2 = double(k+ksv+ksk+ksu), t1_5 = double(n.ids)) 
  # cat.print(sum(tymch$lnls))
  # return(tymch$grad)
  return(matrix(tymch$grad, byrow = FALSE, ncol = length(theta)))
}

.gr.sn.exp.by.i.R <- function(theta, myprod, y, x, zsv, zsk, zsu, n.ids, k, ksv, ksk, ksu, threads = 1, ...) {
  
  eps0  <- y - x %*% theta[1 :k]
  sv    <- sqrt(exp( zsv %*% theta[(k + 1) :(k + ksv)]))
  alp   <- zsk %*% theta[(k + ksv + 1) :(k + ksv + ksk)]
  lam   <- 1 / sqrt(exp( zsu %*% theta[(k + ksv + ksk  + 1) :(k + ksv + ksk + ksu)] ))
  su    <- 1/ lam
  
  epsr  <- eps0 + sv * sqrt(2/pi) * alp / sqrt(1 + alp^2)
  u1    <- (myprod*epsr + lam*sv^2)/sv
  b     <- alp * myprod
  a     <- -b * lam * sv
  a2    <- a / sqrt(1 + b^2)
  
  term1 <- -TOwen( u1, a2 / u1, threads = 1 )
  term2 <- -TOwen( a2, u1 / a2, threads = 1 )
  term3 <- TOwen( u1, b + a / u1, threads = 1 )
  term4 <- TOwen( a2, b + u1 * (1 + b^2) / a, threads = 1 )
  term5 <- pnorm(a2) * pnorm(-u1)
  t1_5  <- term1 + term2 + term3 + term4 + term5
  
  eb    <- -x
  SVU   <- sv * lam
  SVUv  <- sweep(zsv, 1, SVU * 0.5, FUN = "*") #  SVU * 0.5 * zsv
  SVUu  <- sweep(zsu, 1, -SVU * 0.5, FUN = "*") # -SVU * 0.5 * zsu
  A     <- alp / sqrt(1+alp^2)
  As    <- sweep(zsk, 1, (1+alp^2)^(-3/2), FUN = "*") #zsk * (1-A^2) / sqrt(1+alp^2)
  
  d_u1_b  <- u1b <- sweep(eb, 1, myprod / sv, FUN = "*") # myprod/sv*eb
  d_u1_v  <- u1v <- sweep(zsv, 1, myprod * eps0 * -0.5 / sv, FUN = "*") + SVUv # myprod * eps0 * -0.5 * sv * zsv
  d_u1_s  <- u1s <- myprod * sqrt(2/pi) * As
  d_u1_u  <- u1u <- SVUu
  
  d_a2_b  <- matrix(0, nrow = nrow(x), ncol = ncol(x))
  d_a2_v  <- a2v <- sweep(SVUv, 1, -myprod * A, FUN = "*")  
  d_a2_s  <- a2s <- sweep(As,   1, -myprod * SVU, FUN = "*")  
  d_a2_u  <- a2u <- sweep(SVUu, 1, -myprod * A, FUN = "*")  
  
  d_a2u1_b <-                                  - sweep(u1b, 1, a2 / u1^2, FUN = "*")
  d_a2u1_v <- sweep(a2v, 1, 1 / u1, FUN = "*") - sweep(u1v, 1, a2 / u1^2, FUN = "*")
  d_a2u1_s <- sweep(a2s, 1, 1 / u1, FUN = "*") - sweep(u1s, 1, a2 / u1^2, FUN = "*")
  d_a2u1_u <- sweep(a2u, 1, 1 / u1, FUN = "*") - sweep(u1u, 1, a2 / u1^2, FUN = "*")
  
  d_u1a2_b <- sweep(u1b, 1, 1 / a2, FUN = "*")
  d_u1a2_v <- sweep(u1v, 1, 1 / a2, FUN = "*") - sweep(a2v, 1, u1 / a2^2, FUN = "*")
  d_u1a2_s <- sweep(u1s, 1, 1 / a2, FUN = "*") - sweep(a2s, 1, u1 / a2^2, FUN = "*")
  d_u1a2_u <- sweep(u1u, 1, 1 / a2, FUN = "*") - sweep(a2u, 1, u1 / a2^2, FUN = "*")
  
  av       <- sweep(SVUv, 1, -myprod * alp, FUN = "*")
  alps     <- zsk
  as       <- sweep(alps, 1, -myprod * SVU, FUN = "*")
  au       <- sweep(SVUu, 1, -myprod * alp, FUN = "*")
  
  d_bau1_b <-                                 - sweep(u1b, 1, a / u1^2, FUN = "*")
  d_bau1_v <- sweep(av, 1, 1 / u1, FUN = "*") - sweep(u1v, 1, a / u1^2, FUN = "*")
  d_bau1_s <- sweep(as, 1, 1 / u1, FUN = "*") - sweep(u1s, 1, a / u1^2, FUN = "*") + myprod * alps
  d_bau1_u <- sweep(au, 1, 1 / u1, FUN = "*") - sweep(u1u, 1, a / u1^2, FUN = "*")
  
  d_bu1a_b <- sweep(u1b, 1, (1+b^2) / a, FUN = "*")
  d_bu1a_v <- sweep(u1v, 1, (1+b^2) / a, FUN = "*") - sweep(av, 1, u1*(1+b^2) / a^2, FUN = "*")
  d_bu1a_u <- sweep(u1u, 1, (1+b^2) / a, FUN = "*") - sweep(au, 1, u1*(1+b^2) / a^2, FUN = "*")
  # d_bu1a_s <- sweep(u1s, 1, (1+b^2) / a, FUN = "*") + sweep(alps, 1, u1 * ( a*2*alp + myprod * SVU *(1+b^2) ) / a^2, FUN = "*") + myprod * alps
  # d_bu1a_s <- sweep(u1s, 1, (1+b^2) / a, FUN = "*") + sweep(alps, 1, myprod*2*b*u1/a, FUN = "*") - sweep(as, 1, u1*(1+b^2)/a^2, FUN = "*") + myprod * alps
  # d_bu1a_s <- sweep(u1s, 1, (1+b^2) / a, FUN = "*") - sweep(as, 1, u1*((1+b^2)/a^2 - 2/SVU^2), FUN = "*") + myprod * alps
  # d_bu1a_s <- sweep(u1s, 1, (1+b^2) / a, FUN = "*") - sweep(as, 1, u1*((1+b^2)/a^2 + 2*b/a/SVU), FUN = "*") + myprod * alps
  d_bu1a_s <- sweep(u1s, 1, (1+b^2) / a, FUN = "*") - sweep(as, 1, u1/SVU^2*(1/b^2-1), FUN = "*") + myprod * alps
  
  # d_ln2su_u <- -0.5*zsu
  
  d_rest_b <- sweep(eb,   1, myprod/su, FUN = "*")
  d_rest_v <- sweep(SVUv, 1, myprod*sqrt(2/pi)*A + SVU, FUN = "*")
  d_rest_s <- sweep(As,   1, myprod*sqrt(2/pi)*SVU, FUN = "*")
  d_rest_u <- sweep(SVUu, 1, myprod*sqrt(2/pi)*A + SVU, FUN = "*") + sweep(zsu, 1, myprod * eps0 * -0.5 / su -0.5, FUN = "*") #- 0.5*zsu
  
  T_gr_1_u1_a2u1 <- TOwen.gr.1(u1,a2/u1); T_gr_2_u1_a2u1 <- TOwen.gr.2(u1,a2/u1)
  
  T_gr_1_a2_u1a2 <- TOwen.gr.1(a2,u1/a2); T_gr_2_a2_u1a2 <- TOwen.gr.2(a2,u1/a2)
  
  T_gr_1_u1_bau1 <- TOwen.gr.1(u1,b+a/u1); T_gr_2_u1_bau1 <- TOwen.gr.2(u1,b+a/u1)
  
  T_gr_1_a2_bu1a <- TOwen.gr.1(a2,b+u1*(1+b^2)/a); T_gr_2_a2_bu1a <- TOwen.gr.2(a2,b+u1*(1+b^2)/a)
  
  d_term1_b <- sweep(d_u1_b, 1, T_gr_1_u1_a2u1, FUN = "*") + sweep(d_a2u1_b, 1, T_gr_2_u1_a2u1, FUN = "*") 
  d_term2_b <- sweep(d_a2_b, 1, T_gr_1_a2_u1a2, FUN = "*") + sweep(d_u1a2_b, 1, T_gr_2_a2_u1a2, FUN = "*") 
  d_term3_b <- sweep(d_u1_b, 1, T_gr_1_u1_bau1, FUN = "*") + sweep(d_bau1_b, 1, T_gr_2_u1_bau1, FUN = "*") 
  d_term4_b <- sweep(d_a2_b, 1, T_gr_1_a2_bu1a, FUN = "*") + sweep(d_bu1a_b, 1, T_gr_2_a2_bu1a, FUN = "*")
  
  d_term1_v <- sweep(d_u1_v, 1, T_gr_1_u1_a2u1, FUN = "*") + sweep(d_a2u1_v, 1, T_gr_2_u1_a2u1, FUN = "*") 
  d_term2_v <- sweep(d_a2_v, 1, T_gr_1_a2_u1a2, FUN = "*") + sweep(d_u1a2_v, 1, T_gr_2_a2_u1a2, FUN = "*") 
  d_term3_v <- sweep(d_u1_v, 1, T_gr_1_u1_bau1, FUN = "*") + sweep(d_bau1_v, 1, T_gr_2_u1_bau1, FUN = "*") 
  d_term4_v <- sweep(d_a2_v, 1, T_gr_1_a2_bu1a, FUN = "*") + sweep(d_bu1a_v, 1, T_gr_2_a2_bu1a, FUN = "*")
  
  d_term1_s <- sweep(d_u1_s, 1, T_gr_1_u1_a2u1, FUN = "*") + sweep(d_a2u1_s, 1, T_gr_2_u1_a2u1, FUN = "*") 
  d_term2_s <- sweep(d_a2_s, 1, T_gr_1_a2_u1a2, FUN = "*") + sweep(d_u1a2_s, 1, T_gr_2_a2_u1a2, FUN = "*") 
  d_term3_s <- sweep(d_u1_s, 1, T_gr_1_u1_bau1, FUN = "*") + sweep(d_bau1_s, 1, T_gr_2_u1_bau1, FUN = "*") 
  d_term4_s <- sweep(d_a2_s, 1, T_gr_1_a2_bu1a, FUN = "*") + sweep(d_bu1a_s, 1, T_gr_2_a2_bu1a, FUN = "*")
  
  d_term1_u <- sweep(d_u1_u, 1, T_gr_1_u1_a2u1, FUN = "*") + sweep(d_a2u1_u, 1, T_gr_2_u1_a2u1, FUN = "*") 
  d_term2_u <- sweep(d_a2_u, 1, T_gr_1_a2_u1a2, FUN = "*") + sweep(d_u1a2_u, 1, T_gr_2_a2_u1a2, FUN = "*") 
  d_term3_u <- sweep(d_u1_u, 1, T_gr_1_u1_bau1, FUN = "*") + sweep(d_bau1_u, 1, T_gr_2_u1_bau1, FUN = "*") 
  d_term4_u <- sweep(d_a2_u, 1, T_gr_1_a2_bu1a, FUN = "*") + sweep(d_bu1a_u, 1, T_gr_2_a2_bu1a, FUN = "*")
  
  Phis_gr_1 <- pnorm(a2) * dnorm(-u1)
  Phis_gr_2 <- pnorm(-u1) * dnorm(a2)
  
  d_term5_b <- sweep(-d_u1_b, 1, Phis_gr_1, FUN = "*") + sweep(d_a2_b, 1, Phis_gr_2, FUN = "*")
  d_term5_v <- sweep(-d_u1_v, 1, Phis_gr_1, FUN = "*") + sweep(d_a2_v, 1, Phis_gr_2, FUN = "*")
  d_term5_s <- sweep(-d_u1_s, 1, Phis_gr_1, FUN = "*") + sweep(d_a2_s, 1, Phis_gr_2, FUN = "*")
  d_term5_u <- sweep(-d_u1_u, 1, Phis_gr_1, FUN = "*") + sweep(d_a2_u, 1, Phis_gr_2, FUN = "*")
  
  g_b <- d_rest_b + sweep(-d_term1_b-d_term2_b+d_term3_b+d_term4_b+d_term5_b, 1, t1_5, FUN = "/")
  g_v <- d_rest_v + sweep(-d_term1_v-d_term2_v+d_term3_v+d_term4_v+d_term5_v, 1, t1_5, FUN = "/")
  g_s <- d_rest_s + sweep(-d_term1_s-d_term2_s+d_term3_s+d_term4_s+d_term5_s, 1, t1_5, FUN = "/")
  g_u <- d_rest_u + sweep(-d_term1_u-d_term2_u+d_term3_u+d_term4_u+d_term5_u, 1, t1_5, FUN = "/")
  
  g0  <- cbind(g_b, g_v, g_s, g_u)
  colnames(g0) <- names(theta)
  return(ifelse(is.na(g0) | is.infinite(g0), 0, g0))
  
}


# hessian -----------------------------------------------------------------

.hess.sn.exp <- function(theta, myprod, y, x, zsv, zsk, zsu, n.ids, k, ksv, ksk, ksu, threads = 1, ...) {
  # cat.print(theta)  
  h0 <- tryCatch(
    hessian2(funcg = .gr.sn.exp, at = theta, myprod=myprod, y=y, x=x, zsv=zsv, zsk=zsk, zsu=zsu, n.ids=n.ids, k=k, ksv=ksv, ksk=ksk, ksu=ksu, threads = threads),
    error = function(e) e )
  if(!inherits(h0, "error")){
    return(ifelse(is.na(h0) | is.infinite(h0), 0, h0) )
  } else {
    return(NULL)
  }
  
}

.hess.sn.exp.bhhh <- function(theta, myprod, y, x, zsv, zsk, zsu, n.ids, k, ksv, ksk, ksu, threads = 1, ...) {
  # cat("this is .hess.sn.exp.bhhh \n\n")
  # cat.print(theta)
  g0 <- .gr.sn.exp.by.i(theta, myprod=myprod, y=y, x=x, zsv=zsv, zsk=zsk, zsu=zsu, n.ids=n.ids, k=k, ksv=ksv, ksk=ksk, ksu=ksu, threads = threads)
  # cat.print(g0)
  h0 <- -crossprod(g0)
  # cat.print(h0)
  return(ifelse(is.na(h0) | is.infinite(h0), 0, h0) )
  
}


# grad + hess -------------------------------------------------------------

.gr.hess.sn.exp <- function(theta, myprod, y, x, zsv, zsk, zsu, n.ids, k, ksv, ksk, ksu, threads = 1, ...) {
  tymch <- .C("sn_exp_ll_grad", as.integer(threads),
              as.double(myprod), as.double(y), as.double(x),
              as.double(zsv), as.double(zsk), as.double(zsu),
              as.integer(n.ids), as.integer(k), 
              as.integer(ksv), as.integer(ksk), as.integer(ksu),
              as.double(theta), 
              lnls = double(n.ids), lnl = as.double(1), 
              grad = double(n.ids*(k+ksv+ksk+ksu)), grad2 = double(k+ksv+ksk+ksu), t1_5 = double(n.ids)) 
  
  h0 <- tryCatch(
    hessian2(funcg = .gr.sn.exp, at = theta, myprod=myprod, y=y, x=x, zsv=zsv, zsk=zsk, zsu=zsu, n.ids=n.ids, k=k, ksv=ksv, ksk=ksk, ksu=ksu, threads = threads),
    error = function(e) e )
  
  list( grad = ifelse(is.na(tymch$grad2) | is.infinite(tymch$grad2), 0, tymch$grad2),  hessian1 = ifelse(is.na(h0) | is.infinite(h0), 0, h0) )
  
  
}

.gr.hess.sn.exp_ <- function(theta, myprod, y, x, zsv, zsk, zsu, n.ids, k, ksv, ksk, ksu, threads = 1, ...) {
  
  # g0 <- gradient1(func = .ll.sn.exp, at = theta, myprod=myprod, y=y, x=x, zsv=zsv, zsk=zsk, zsu=zsu, n.ids=n.ids, k=k, ksv=ksv, ksk=ksk, ksu=ksu, threads = threads)
  eps0  <- y - x %*% theta[1 :k]
  sv    <- sqrt(exp( zsv %*% theta[(k + 1) :(k + ksv)]))
  alp   <- zsk %*% theta[(k + ksv + 1) :(k + ksv + ksk)]
  lam   <- 1 / sqrt(exp( zsu %*% theta[(k + ksv + ksk  + 1) :(k + ksv + ksk + ksu)] ))
  su    <- 1/ lam
  
  epsr  <- eps0 + sv * sqrt(2/pi) * alp / sqrt(1 + alp^2)
  u1    <- (myprod*epsr + lam*sv^2)/sv
  b     <- alp * myprod
  a     <- -b * lam * sv
  a2    <- a / sqrt(1 + b^2)
  
  term1 <- -TOwen( u1, a2 / u1, threads = 1 )
  term2 <- -TOwen( a2, u1 / a2, threads = 1 )
  term3 <- TOwen( u1, b + a / u1, threads = 1 )
  term4 <- TOwen( a2, b + u1 * (1 + b^2) / a, threads = 1 )
  term5 <- pnorm(a2) * pnorm(-u1)
  t1_5  <- term1 + term2 + term3 + term4 + term5
  
  eb    <- -x
  SVU   <- sv * lam
  SVUv  <- sweep(zsv, 1, SVU * 0.5, FUN = "*") #  SVU * 0.5 * zsv
  SVUu  <- sweep(zsu, 1, -SVU * 0.5, FUN = "*") # -SVU * 0.5 * zsu
  A     <- alp / sqrt(1+alp^2)
  As    <- sweep(zsk, 1, (1+alp^2)^(-3/2), FUN = "*") #zsk * (1-A^2) / sqrt(1+alp^2)
  
  d_u1_b  <- u1b <- sweep(eb, 1, myprod / sv, FUN = "*") # myprod/sv*eb
  d_u1_v  <- u1v <- sweep(zsv, 1, myprod * eps0 * -0.5 / sv, FUN = "*") + SVUv # myprod * eps0 * -0.5 * sv * zsv
  d_u1_s  <- u1s <- myprod * sqrt(2/pi) * As
  d_u1_u  <- u1u <- SVUu
  
  d_a2_b  <- matrix(0, nrow = nrow(x), ncol = ncol(x))
  d_a2_v  <- a2v <- sweep(SVUv, 1, -myprod * A, FUN = "*")  
  d_a2_s  <- a2s <- sweep(As,   1, -myprod * SVU, FUN = "*")  
  d_a2_u  <- a2u <- sweep(SVUu, 1, -myprod * A, FUN = "*")  
  
  d_a2u1_b <-                                  - sweep(u1b, 1, a2 / u1^2, FUN = "*")
  d_a2u1_v <- sweep(a2v, 1, 1 / u1, FUN = "*") - sweep(u1v, 1, a2 / u1^2, FUN = "*")
  d_a2u1_s <- sweep(a2s, 1, 1 / u1, FUN = "*") - sweep(u1s, 1, a2 / u1^2, FUN = "*")
  d_a2u1_u <- sweep(a2u, 1, 1 / u1, FUN = "*") - sweep(u1u, 1, a2 / u1^2, FUN = "*")
  
  d_u1a2_b <- sweep(u1b, 1, 1 / a2, FUN = "*")
  d_u1a2_v <- sweep(u1v, 1, 1 / a2, FUN = "*") - sweep(a2v, 1, u1 / a2^2, FUN = "*")
  d_u1a2_s <- sweep(u1s, 1, 1 / a2, FUN = "*") - sweep(a2s, 1, u1 / a2^2, FUN = "*")
  d_u1a2_u <- sweep(u1u, 1, 1 / a2, FUN = "*") - sweep(a2u, 1, u1 / a2^2, FUN = "*")
  
  av       <- sweep(SVUv, 1, -myprod * alp, FUN = "*")
  alps     <- zsk
  as       <- sweep(alps, 1, -myprod * SVU, FUN = "*")
  au       <- sweep(SVUu, 1, -myprod * alp, FUN = "*")
  
  d_bau1_b <-                                 - sweep(u1b, 1, a / u1^2, FUN = "*")
  d_bau1_v <- sweep(av, 1, 1 / u1, FUN = "*") - sweep(u1v, 1, a / u1^2, FUN = "*")
  d_bau1_s <- sweep(as, 1, 1 / u1, FUN = "*") - sweep(u1s, 1, a / u1^2, FUN = "*") + myprod * alps
  d_bau1_u <- sweep(au, 1, 1 / u1, FUN = "*") - sweep(u1u, 1, a / u1^2, FUN = "*")
  
  d_bu1a_b <- sweep(u1b, 1, (1+b^2) / a, FUN = "*")
  d_bu1a_v <- sweep(u1v, 1, (1+b^2) / a, FUN = "*") - sweep(av, 1, u1*(1+b^2) / a^2, FUN = "*")
  d_bu1a_u <- sweep(u1u, 1, (1+b^2) / a, FUN = "*") - sweep(au, 1, u1*(1+b^2) / a^2, FUN = "*")
  # d_bu1a_s <- sweep(u1s, 1, (1+b^2) / a, FUN = "*") + sweep(alps, 1, u1 * ( a*2*alp + myprod * SVU *(1+b^2) ) / a^2, FUN = "*") + myprod * alps
  # d_bu1a_s <- sweep(u1s, 1, (1+b^2) / a, FUN = "*") + sweep(alps, 1, myprod*2*b*u1/a, FUN = "*") - sweep(as, 1, u1*(1+b^2)/a^2, FUN = "*") + myprod * alps
  # d_bu1a_s <- sweep(u1s, 1, (1+b^2) / a, FUN = "*") - sweep(as, 1, u1*((1+b^2)/a^2 - 2/SVU^2), FUN = "*") + myprod * alps
  # d_bu1a_s <- sweep(u1s, 1, (1+b^2) / a, FUN = "*") - sweep(as, 1, u1*((1+b^2)/a^2 + 2*b/a/SVU), FUN = "*") + myprod * alps
  d_bu1a_s <- sweep(u1s, 1, (1+b^2) / a, FUN = "*") - sweep(as, 1, u1/SVU^2*(1/b^2-1), FUN = "*") + myprod * alps
  
  # d_ln2su_u <- -0.5*zsu
  
  d_rest_b <- sweep(eb,   1, myprod/su, FUN = "*")
  d_rest_v <- sweep(SVUv, 1, myprod*sqrt(2/pi)*A + SVU, FUN = "*")
  d_rest_s <- sweep(As,   1, myprod*sqrt(2/pi)*SVU, FUN = "*")
  d_rest_u <- sweep(SVUu, 1, myprod*sqrt(2/pi)*A + SVU, FUN = "*") + sweep(zsu, 1, myprod * eps0 * -0.5 / su -0.5, FUN = "*") #- 0.5*zsu
  
  T_gr_1_u1_a2u1 <- TOwen.gr.1(u1,a2/u1); T_gr_2_u1_a2u1 <- TOwen.gr.2(u1,a2/u1)
  
  T_gr_1_a2_u1a2 <- TOwen.gr.1(a2,u1/a2); T_gr_2_a2_u1a2 <- TOwen.gr.2(a2,u1/a2)
  
  T_gr_1_u1_bau1 <- TOwen.gr.1(u1,b+a/u1); T_gr_2_u1_bau1 <- TOwen.gr.2(u1,b+a/u1)
  
  T_gr_1_a2_bu1a <- TOwen.gr.1(a2,b+u1*(1+b^2)/a); T_gr_2_a2_bu1a <- TOwen.gr.2(a2,b+u1*(1+b^2)/a)
  
  d_term1_b <- sweep(d_u1_b, 1, T_gr_1_u1_a2u1, FUN = "*") + sweep(d_a2u1_b, 1, T_gr_2_u1_a2u1, FUN = "*") 
  d_term2_b <- sweep(d_a2_b, 1, T_gr_1_a2_u1a2, FUN = "*") + sweep(d_u1a2_b, 1, T_gr_2_a2_u1a2, FUN = "*") 
  d_term3_b <- sweep(d_u1_b, 1, T_gr_1_u1_bau1, FUN = "*") + sweep(d_bau1_b, 1, T_gr_2_u1_bau1, FUN = "*") 
  d_term4_b <- sweep(d_a2_b, 1, T_gr_1_a2_bu1a, FUN = "*") + sweep(d_bu1a_b, 1, T_gr_2_a2_bu1a, FUN = "*")
  
  d_term1_v <- sweep(d_u1_v, 1, T_gr_1_u1_a2u1, FUN = "*") + sweep(d_a2u1_v, 1, T_gr_2_u1_a2u1, FUN = "*") 
  d_term2_v <- sweep(d_a2_v, 1, T_gr_1_a2_u1a2, FUN = "*") + sweep(d_u1a2_v, 1, T_gr_2_a2_u1a2, FUN = "*") 
  d_term3_v <- sweep(d_u1_v, 1, T_gr_1_u1_bau1, FUN = "*") + sweep(d_bau1_v, 1, T_gr_2_u1_bau1, FUN = "*") 
  d_term4_v <- sweep(d_a2_v, 1, T_gr_1_a2_bu1a, FUN = "*") + sweep(d_bu1a_v, 1, T_gr_2_a2_bu1a, FUN = "*")
  
  d_term1_s <- sweep(d_u1_s, 1, T_gr_1_u1_a2u1, FUN = "*") + sweep(d_a2u1_s, 1, T_gr_2_u1_a2u1, FUN = "*") 
  d_term2_s <- sweep(d_a2_s, 1, T_gr_1_a2_u1a2, FUN = "*") + sweep(d_u1a2_s, 1, T_gr_2_a2_u1a2, FUN = "*") 
  d_term3_s <- sweep(d_u1_s, 1, T_gr_1_u1_bau1, FUN = "*") + sweep(d_bau1_s, 1, T_gr_2_u1_bau1, FUN = "*") 
  d_term4_s <- sweep(d_a2_s, 1, T_gr_1_a2_bu1a, FUN = "*") + sweep(d_bu1a_s, 1, T_gr_2_a2_bu1a, FUN = "*")
  
  d_term1_u <- sweep(d_u1_u, 1, T_gr_1_u1_a2u1, FUN = "*") + sweep(d_a2u1_u, 1, T_gr_2_u1_a2u1, FUN = "*") 
  d_term2_u <- sweep(d_a2_u, 1, T_gr_1_a2_u1a2, FUN = "*") + sweep(d_u1a2_u, 1, T_gr_2_a2_u1a2, FUN = "*") 
  d_term3_u <- sweep(d_u1_u, 1, T_gr_1_u1_bau1, FUN = "*") + sweep(d_bau1_u, 1, T_gr_2_u1_bau1, FUN = "*") 
  d_term4_u <- sweep(d_a2_u, 1, T_gr_1_a2_bu1a, FUN = "*") + sweep(d_bu1a_u, 1, T_gr_2_a2_bu1a, FUN = "*")
  
  Phis_gr_1 <- pnorm(a2) * dnorm(-u1)
  Phis_gr_2 <- pnorm(-u1) * dnorm(a2)
  
  d_term5_b <- sweep(-d_u1_b, 1, Phis_gr_1, FUN = "*") + sweep(d_a2_b, 1, Phis_gr_2, FUN = "*")
  d_term5_v <- sweep(-d_u1_v, 1, Phis_gr_1, FUN = "*") + sweep(d_a2_v, 1, Phis_gr_2, FUN = "*")
  d_term5_s <- sweep(-d_u1_s, 1, Phis_gr_1, FUN = "*") + sweep(d_a2_s, 1, Phis_gr_2, FUN = "*")
  d_term5_u <- sweep(-d_u1_u, 1, Phis_gr_1, FUN = "*") + sweep(d_a2_u, 1, Phis_gr_2, FUN = "*")
  
  g_b <- d_rest_b + sweep(-d_term1_b-d_term2_b+d_term3_b+d_term4_b+d_term5_b, 1, t1_5, FUN = "/")
  g_v <- d_rest_v + sweep(-d_term1_v-d_term2_v+d_term3_v+d_term4_v+d_term5_v, 1, t1_5, FUN = "/")
  g_s <- d_rest_s + sweep(-d_term1_s-d_term2_s+d_term3_s+d_term4_s+d_term5_s, 1, t1_5, FUN = "/")
  g_u <- d_rest_u + sweep(-d_term1_u-d_term2_u+d_term3_u+d_term4_u+d_term5_u, 1, t1_5, FUN = "/")
  
  g0  <- colSums(cbind(g_b, g_v, g_s, g_u), na.rm = TRUE)
  
  names(g0) <- names(theta)
  h0 <- .hess.sn.exp.bhhh(theta, myprod=myprod, y=y, x=x, zsv=zsv, zsk=zsk, zsu=zsu, n.ids=n.ids, k=k, ksv=ksv, ksk=ksk, ksu=ksu, threads = threads)
  
  list( grad = ifelse(is.na(g0) | is.infinite(g0), 0, g0),  hessian1 = ifelse(is.na(h0) | is.infinite(h0), 0, h0) )
  
}

.gr.hess.sn.exp.bhhh <- function(theta, myprod, y, x, zsv, zsk, zsu, n.ids, k, ksv, ksk, ksu, threads = 1, ...) {
  
  g0 <- gradient1(func = .ll.sn.exp, at = theta, myprod=myprod, y=y, x=x, zsv=zsv, zsk=zsk, zsu=zsu, n.ids=n.ids, k=k, ksv=ksv, ksk=ksk, ksu=ksu, threads = threads)
  names(g0) <- names(theta)
  h0 <- hessian1(func = .ll.sn.exp, at = theta, myprod=myprod, y=y, x=x, zsv=zsv, zsk=zsk, zsu=zsu, n.ids=n.ids, k=k, ksv=ksv, ksk=ksk, ksu=ksu, threads = threads)
  
  list( grad = ifelse(is.na(g0) | is.infinite(g0), 0, g0),  hessian1 = ifelse(is.na(h0) | is.infinite(h0), 0, h0) )
  
}


# SN-TN -------------------------------------------------------------------

# log likelihood ----------------------------------------------------------

.ll.sn.tn <- function(theta, myprod, y, x, zsv, zsk, zsu, zmu = NULL, 
                      n.ids, k, ksv, ksk, ksu, kmu = NULL, tn = 1, threads = 1, ...){
  # cat.print(tn)
  tymch <- .C("ll_sn_tn", as.integer(threads),
              as.double(myprod), as.double(y), as.double(x),
              as.double(zsv), as.double(zsk), as.double(zsu), as.double(zmu),
              as.integer(n.ids), as.integer(k), 
              as.integer(ksv), as.integer(ksk), as.integer(ksu), as.integer(kmu),
              as.double(theta), as.double(tn),
              lnls = double(n.ids), lnl = as.double(1)) 
  # cat.print(sum(tymch$lnls))
  return(tymch$lnl)
}

.ll.sn.tn.R <- function(theta, myprod, y, x, zsv, zsk, zsu, zmu = matrix(0,n.ids,1), 
                        n.ids, k, ksv, ksk, ksu, kmu = 0, threads = 1, ...){
  # k:    betas
  # ksv:  noise, scale/variance
  # ksk:  noise, skew
  # ksu:  inefficiency: scale/variance
  # kmu:  inefficiency: location
  # cat.print(1)
  
  eps0  <- y - x %*%    theta[1                        :k]
  # cat.print(2)
  sv2   <- exp( zsv %*% theta[(k + 1)                  :(k + ksv)])
  # cat.print(3)
  al0   <- zsk %*%      theta[(k + ksv + 1)            :(k + ksv + ksk)]
  # cat.print(4)
  su2   <- exp( zsu %*% theta[(k + ksv + ksk  + 1)     :(k + ksv + ksk + ksu)] )
  # cat.print(5)
  mu0   <- zmu %*%      theta[(k + ksv + ksk + ksu + 1):(k + ksv + ksk + ksu + kmu)]
  # cat.print(6)
  
  sv    <- sqrt(sv2)
  su    <- sqrt(su2)
  sig   <- sqrt(sv2 + su2)
  sstar <- sv*su/sig
  epsr  <- eps0 + sv * sqrt(2/pi) * al0 / sqrt(1+al0^2)
  mu1   <- (mu0*sv2 - epsr*myprod*su2) / sig^2
  b1    <- al0 / sv * myprod * sstar 
  a1    <- al0 / sv * (epsr + myprod * mu1)
  a2    <- a1 / sqrt(1+b1^2)
  u1    <- -mu1 / sstar
  
  # if(any(is.na(a2/u1))) print(as.vector(a2/u1))
  # 
  # if(sum(is.na(a2/u1)) > 0){
  #   tymch1 <- data.frame(eps0=unname(eps0),sv,al0,su,mu0,epsr=unname(epsr),a1,b1,mu1=unname(mupr),a2=unname(a2))
  #   cat.print(head(tymch1,17))
  #   cat.print(theta)
  #   warning("Something is wrong, NA")
  #   return(tymch1)
  # }
  
  term1 <- -TOwen( u1, a2 / u1, threads )
  term2 <- -TOwen( a2, u1 / a2, threads )
  term3 <- TOwen( u1, b1 + a1 / u1, threads )
  term4 <- TOwen( a2, b1 + u1 * (1 + b1^2) / a1, threads )
  term5 <- pnorm(a2) * pnorm(-u1)
  t1_5  <- term1 + term2 + term3 + term4 + term5
  
  t1_5  <- ifelse(t1_5 < 0, .Machine$double.xmin, t1_5)
  
  tymch1 <- log( t1_5 ) + log(2) - pnorm( mu0/su, log.p = TRUE) - log(sig) +
    dnorm((mu0 + myprod*epsr)/sig, log = TRUE)
  # cat.print(as.vector(tymch1))
  return( sum(tymch1, na.rm = TRUE) )
}

# gradient ----------------------------------------------------------------

.gr.sn.tn <- function(theta, myprod, y, x, zsv, zsk, zsu, zmu = NULL, 
                      n.ids, k, ksv, ksk, ksu, kmu = NULL, tn, threads = 1, ...){
  g0 <- gradient1(func = .ll.sn.tn, at = theta, myprod=myprod, y=y, x=x, zsv=zsv, zsk=zsk, zsu=zsu, zmu = zmu, n.ids=n.ids, k=k, ksv=ksv, ksk=ksk, ksu=ksu, kmu = kmu, tn = tn, threads = threads)
  names(g0) <- names(theta)
  return(ifelse(is.na(g0) | is.infinite(g0), 0, g0))
  
}


# hessian -----------------------------------------------------------------

.hess.sn.tn <- function(theta, myprod, y, x, zsv, zsk, zsu, zmu = NULL, 
                        n.ids, k, ksv, ksk, ksu, kmu = NULL, tn, threads = 1, ...){
  
  h0 <- hessian1(func = .ll.sn.tn, at = theta, myprod=myprod, y=y, x=x, zsv=zsv, zsk=zsk, zsu=zsu, zmu = zmu, n.ids=n.ids, k=k, ksv=ksv, ksk=ksk, ksu=ksu, kmu = kmu, tn = tn, threads = threads)
  return(ifelse(is.na(h0) | is.infinite(h0), 0, h0) )
  
}


# grad + hess -------------------------------------------------------------

.gr.hess.sn.tn <- function(theta, myprod, y, x, zsv, zsk, zsu, zmu = NULL, 
                           n.ids, k, ksv, ksk, ksu, kmu = NULL, tn, threads = 1, ...){
  
  g0 <- gradient1(func = .ll.sn.tn, at = theta, myprod=myprod, y=y, x=x, zsv=zsv, zsk=zsk, zsu=zsu, zmu = zmu, n.ids=n.ids, k=k, ksv=ksv, ksk=ksk, ksu=ksu, kmu = kmu, tn = tn, threads = threads)
  names(g0) <- names(theta)
  h0 <- hessian1(func = .ll.sn.tn, at = theta, myprod=myprod, y=y, x=x, zsv=zsv, zsk=zsk, zsu=zsu, zmu = zmu, n.ids=n.ids, k=k, ksv=ksv, ksk=ksk, ksu=ksu, kmu = kmu, tn = tn, threads = threads)
  
  list( grad = ifelse(is.na(g0) | is.infinite(g0), 0, g0),  hessian1 = ifelse(is.na(h0) | is.infinite(h0), 0, h0) )
  
}

# SN-hN -------------------------------------------------------------------

# log likelihood ----------------------------------------------------------

.ll.sn.hn <- function(theta, myprod, y, x, zsv, zsk, zsu, 
                      n.ids, k, ksv, ksk, ksu, threads = 1, ...){
  # k:    betas
  # ksv:  noise, scale/variance
  # ksk:  noise, skew
  # ksu:  inefficiency: scale/variance
  # kmu:  inefficiency: location
  # cat.print(1)
  
  eps0  <- y - x %*%    theta[1                        :k]
  # cat.print(2)
  sv2   <- exp( zsv %*% theta[(k + 1)                  :(k + ksv)])
  # cat.print(3)
  al0   <- zsk %*%      theta[(k + ksv + 1)            :(k + ksv + ksk)]
  # cat.print(4)
  su2   <- exp( zsu %*% theta[(k + ksv + ksk  + 1)     :(k + ksv + ksk + ksu)] )
  # cat.print(5)
  mu0   <- 0#zmu %*%      theta[(k + ksv + ksk + ksu + 1):(k + ksv + ksk + ksu + kmu)]
  # cat.print(6)
  
  sv    <- sqrt(sv2)
  su    <- sqrt(su2)
  sig   <- sqrt(sv2 + su2)
  sstar <- sv*su/sig
  epsr  <- eps0 + sv * sqrt(2/pi) * al0 / sqrt(1+al0^2)
  mu1   <- (mu0*sv2 - epsr*myprod*su2) / sig^2
  b1    <- al0 / sv * myprod * sstar 
  a1    <- al0 / sv * (epsr + myprod * mu1)
  a2    <- a1 / sqrt(1+b1^2)
  u1    <- -mu1 / sstar
  
  # if(any(is.na(a2/u1))) print(as.vector(a2/u1))
  # 
  # if(sum(is.na(a2/u1)) > 0){
  #   tymch1 <- data.frame(eps0=unname(eps0),sv,al0,su,mu0,epsr=unname(epsr),a1,b1,mu1=unname(mupr),a2=unname(a2))
  #   cat.print(head(tymch1,17))
  #   cat.print(theta)
  #   warning("Something is wrong, NA")
  #   return(tymch1)
  # }
  
  term1 <- -TOwen( u1, a2 / u1, threads )
  term2 <- -TOwen( a2, u1 / a2, threads )
  term3 <- TOwen( u1, b1 + a1 / u1, threads )
  term4 <- TOwen( a2, b1 + u1 * (1 + b1^2) / a1, threads )
  term5 <- pnorm(a2) * pnorm(-u1)
  t1_5  <- term1 + term2 + term3 + term4 + term5
  
  t1_5  <- ifelse(t1_5 < 0, .Machine$double.xmin, t1_5)
  
  tymch1 <- log( t1_5 ) + log(2) - pnorm( mu0/su, log.p = TRUE) - log(sig) +
    dnorm((mu0 + myprod*epsr)/sig, log = TRUE)
  # cat.print(as.vector(tymch1))
  return( sum(tymch1, na.rm = TRUE) )
}

# gradient ----------------------------------------------------------------

.gr.sn.hn <- function(theta, myprod, y, x, zsv, zsk, zsu, n.ids, k, ksv, ksk, ksu, threads = 1, ...) {
  
  g0 <- gradient1(func = .ll.sn.hn, at = theta, myprod=myprod, y=y, x=x, zsv=zsv, zsk=zsk, zsu=zsu, n.ids=n.ids, k=k, ksv=ksv, ksk=ksk, ksu=ksu, threads = threads)
  names(g0) <- names(theta)
  return(ifelse(is.na(g0) | is.infinite(g0), 0, g0))
  
}


# hessian -----------------------------------------------------------------

.hess.sn.hn <- function(theta, myprod, y, x, zsv, zsk, zsu, n.ids, k, ksv, ksk, ksu, threads = 1, ...) {
  
  h0 <- hessian1(func = .ll.sn.hn, at = theta, myprod=myprod, y=y, x=x, zsv=zsv, zsk=zsk, zsu=zsu, n.ids=n.ids, k=k, ksv=ksv, ksk=ksk, ksu=ksu, threads = threads)
  return(ifelse(is.na(h0) | is.infinite(h0), 0, h0) )
  
}


# grad + hess -------------------------------------------------------------

.gr.hess.sn.hn <- function(theta, myprod, y, x, zsv, zsk, zsu, n.ids, k, ksv, ksk, ksu, threads = 1, ...) {
  
  g0 <- gradient1(func = .ll.sn.hn, at = theta, myprod=myprod, y=y, x=x, zsv=zsv, zsk=zsk, zsu=zsu, n.ids=n.ids, k=k, ksv=ksv, ksk=ksk, ksu=ksu, threads = threads)
  names(g0) <- names(theta)
  h0 <- hessian1(func = .ll.sn.hn, at = theta, myprod=myprod, y=y, x=x, zsv=zsv, zsk=zsk, zsu=zsu, n.ids=n.ids, k=k, ksv=ksv, ksk=ksk, ksu=ksu, threads = threads)
  
  list( grad = ifelse(is.na(g0) | is.infinite(g0), 0, g0),  hessian1 = ifelse(is.na(h0) | is.infinite(h0), 0, h0) )
  
}


# Print the estimation results --------------------------------------------

.printoutcs = function(x, digits, k, ksv, ksu, ksk, kmu, na.print = "NA", dist, max.name.length, mycutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), mysymbols = c("***", "**", "*", ".", " ")){
  
  Cf = cbind(ifelse(x[,1, drop = FALSE]> 999, formatC(x[,1, drop = FALSE], digits = 1, format = "e",width = 10), formatC(x[,1, drop = FALSE], digits = digits, format = "f", width = 10)),
             ifelse(x[,2, drop = FALSE]>999, formatC(x[,2, drop = FALSE], digits = 1, format = "e", width = 10), formatC(x[,2, drop = FALSE], digits = digits, format = "f", width = 10)),
             ifelse(x[,3, drop = FALSE]>999, formatC(x[,3, drop = FALSE], digits = 1, format = "e", width = 7), formatC(x[,3, drop = FALSE], digits = 2, format = "f", width = 7)),
             ifelse(x[,4, drop = FALSE]>999, formatC(x[,4, drop = FALSE], digits = 1, format = "e", width = 10), formatC(x[,4, drop = FALSE], digits = digits, format = "f", width = 10)),
             formatC(mysymbols[findInterval(x = x[,4], vec = mycutpoints)], flag = "-"))
  
  # mycutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1)
  # mysymbols = c("***", "**", "*", ".", " ")
  #
  # pvals <- m0$table[,4]
  #
  # findInterval(x = pvals, vec = mycutpoints)
  #
  # cbind(pvals,mysymbols[findInterval(x = pvals, vec = cutpoints)])
  #
  # pval_sym <- mysymbols[findInterval(x = x[,4], vec = cutpoints)]
  
  
  # cat("               Coef.        SE       z       P>|z|\n", sep = "")
  row.names(Cf) <- formatC(row.names(Cf), width = max(nchar(row.names(Cf))), flag = "-")
  cat("",rep(" ", max.name.length+6),"Coef.        SE       z       P>|z|\n", sep = "")
  dimnames(Cf)[[2]] <- rep("", dim(Cf)[[2]])
  cat("",rep("_", max.name.length+42-1),"", "\n", "Frontier", "\n", sep = "")
  print.default(Cf[1:k,,drop=FALSE], quote = FALSE, right = TRUE, na.print = na.print)
  
  cat("",rep("-", max.name.length+42-1),"", "\n", "Variance of the random noise component: log(sigma_v^2)", "\n", sep = "")
  # dimnames(Cf)[[2]] <- rep("", dim(Cf)[[2]])
  print.default(Cf[(k+1):(k+ksv),,drop=FALSE], quote = FALSE, right = TRUE, na.print = na.print)
  
  cat("",rep("-", max.name.length+42-1),"", "\n", "Skewness of the random noise component: `alpha`", "\n", sep = "")
  # dimnames(Cf)[[2]] <- rep("", dim(Cf)[[2]])
  print.default(Cf[(k + ksv + 1):(k + ksv + ksk),,drop=FALSE], quote = FALSE, right = TRUE, na.print = na.print)
  
  cat("",rep("-", max.name.length+42-1),"", "\n", "Inefficiency component: log(sigma_u^2)", "\n", sep = "")
  print.default(Cf[(k + ksv + ksk+1):(k + ksv + ksk+ksu),,drop=FALSE], quote = FALSE, right = TRUE, na.print = na.print)
  if(dist == "t"){
    cat("",rep("-", max.name.length+42-1),"", "\n", "Mu (location parameter)", "\n", sep = "")
    print.default(Cf[(k + ksv + ksk+ksu+1):(k + ksv + ksk+ksu+kmu),,drop=F], quote = FALSE, right = TRUE, na.print = na.print)
  }
  # if(nrow(Cf[-c(1:(k+kv+ku+kdel)),,drop=FALSE]) >= 1){
  #   cat("",rep("-", max.name.length+42-1),"", "\n", "Parameters of compound error distribution", "\n", sep = "")
  #   print.default(Cf[-c(1:(k+kv+ku+kdel)),,drop=FALSE], quote = FALSE, right = TRUE, na.print = na.print)
  # }
  cat("",rep("_", max.name.length+42-1),"", "\n", sep = "")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  invisible(x)
}

# time ----
.timing <- function(x, wording = "Time elapsed is")
{
  if(x < 60){
    # cat("\n")
    cat("",wording,"",x," seconds\n", sep = "")
    # cat("\n")
  } else {
    if(x >= 60 & x < 60*60){
      minutes <- floor(x/60)
      seconds <- round(x - minutes * 60,1)
      # cat("\n")
      cat("",wording,"",minutes," minute(s) and ",seconds," second(s)\n", sep = "")
      # cat("\n")
    } else {
      if(x >= 60*60 & x < 60*60*24){
        hours   <- floor(x / 60 / 60)
        minutes <- round( (x - hours * 60 *60) / 60, 1)
        seconds <- floor(x - hours * 60 *60 - minutes * 60)
        # cat("\n")
        cat("",wording,"",hours," hour(s) and ",minutes," minute(s) \n", sep = "")
        # cat("\n")
      } else {
        if(x >= 60*60*24){
          days    <- floor(x / 60 / 60 / 24)
          hours   <- round( (x - days * 60 * 60 * 24) / 60 /60 ,1)
          minutes <- floor( (x - days * 60 * 60 * 24 - hours * 60 *60) / 60)
          seconds <- floor(x - days * 60 * 60 * 24 - hours * 60 *60 - minutes * 60)
          # cat("\n")
          cat("",wording,"",days," day(s) and ",hours," hour(s)\n", sep = "")
          # cat("\n")
        }
      }
    }
  }
}


.my.prettyNum <- function(xx2){
  n <- length(xx2)
  xx2.pn <- prettyNum(xx2)
  my.integers <- !(1:n %in%  grep(".", xx2.pn, fixed = TRUE))
  xx3 <- ifelse(abs(xx2) < 1e-04 & xx2 != 0, formatC(xx2, digits = 1, format = "e", width = 5), ifelse(abs(xx2) >= 1e-04 & abs(xx2) < 10, formatC(xx2, format = "f", digits = 4), ifelse(abs(xx2) >= 10 & abs(xx2) < 100, formatC(xx2, format = "f", digits = 3), ifelse(abs(xx2) >= 100 & abs(xx2) < 1000, formatC(xx2, format = "f", digits = 2), ifelse(abs(xx2) >= 1000, formatC(xx2, format = "e", digits = 1), formatC(xx2, format = "f", digits = 1))))))
  xx3[my.integers] <- xx2.pn[my.integers]
  xx3
}

# my summary function
.su <- function(x, mat.var.in.col = TRUE, digits = 4, probs = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9), print = FALSE, names = NULL, ...){
  
  xvec2 <- xvec1 <- FALSE
  # print(class(x))
  
  if(is.list(x)){
    # cat("this is a list\n")
    mynames <- paste0("Var",1:length(x))
  }
  
  if(is.matrix(x)){
    # cat("this is a matrix\n")
    if(min(dim(x)) == 1){
      # xvec1 <- TRUE
      if(which(dim(x) == 1) == 2){
        mynames <- colnames(x)
        
      } else {
        mynames <- rownames(x)
        x <- t(x)
      }
      # x <- as.vector(x)
    } else {
      if(!mat.var.in.col){
        x <- t(x)
      }
      mynames <- colnames(x)
    }
    
    # print(mynames)
    if(is.null(mynames)) mynames <- paste("Var", seq_len(ncol(x)), sep = "")
    x <- as.data.frame(x)
    # print(x)
    # mynames <- colnames(x)
  } # end if matrix
  
  if(is.vector(x) & !is.list(x)){
    # cat("this is a vector\n")
    xvec2 <- TRUE
    mynames <- deparse(substitute(x))
    x <- data.frame(Var1 = x)
  } # end if vector
  
  # print(mynames)
  # print(names)
  
  if(!is.null(names)){
    if(length(mynames) == length(names)){
      mynames <- names
    }
  }
  
  # cat("nymanes", sep ="")
  # print(mynames)
  if(is.list(x)){
    t1 <- sapply(x, function(x) c(Obs = length(x), NAs = sum(is.na(x)), Mean = mean(x, na.rm = TRUE), StDev = sd(x, na.rm = TRUE), IQR = IQR(x, na.rm = TRUE), Min = min(x, na.rm = TRUE), quantile( x, probs = probs, na.rm = TRUE ), Max = max(x, na.rm = TRUE)))
    # print(t1)
    # print(mynames)
    # print(class(t1))
    # print(is.integer(t1))
    # print( formatC(t1, format = "f") == formatC(t1, format = "fg") )
    # print(dim(t1))
    # if(xvec2 & !xvec1) colnames(t1) <- mynames
    colnames(t1) <- mynames
    if(print){
      # t2 <- rbind(
      #   formatC(t1[c(1:2),,drop = FALSE], digits = digits, format = "fg"),
      #   formatC(t1[-c(1:2),,drop = FALSE], digits = digits, format = "f")
      # )
      print(data.frame(.my.prettyNum(t1)))
    }
  } else if(!is.vector(x) & !is.matrix(x) & !is.data.frame(x)){
    stop("Provide list, vector, matrix, or data.frame")
  } else {
    t1 <- apply(x, 2, function(x) c(Obs = length(x), NAs = sum(is.na(x)), Mean = mean(x, na.rm = TRUE), StDev = sd(x, na.rm = TRUE), IQR = IQR(x, na.rm = TRUE), Min = min(x, na.rm = TRUE), quantile( x, probs = probs, na.rm = TRUE ), Max = max(x, na.rm = TRUE)))
    # print(t1)
    # print(mynames)
    print(is.integer(t1))
    print(class(t1))
    # print(dim(t1))
    if(xvec2 & !xvec1) colnames(t1) <- mynames
    if(print) print(data.frame(.my.prettyNum(t1)))#t1 <- data.frame(formatC(t1, digits = digits, ...)); print(t1)
  }
  tymch <- t(t1)
  class(tymch) <- "snreg"
  return(tymch)
}

.mlsearch <- function(Formula, logL, init = NULL, int_inds = NULL,
                      auxilary.reports = FALSE, report.rescaling = FALSE,
                      seed = 173451, ...){
  
  set.seed(seed)
  
  # require(Formula)
  formula1 <- Formula
  # print(formula1)
  neqn <- length(formula1)[2]-1
  # print(neqn)
  
  # names for the RHSs
  # determine indicies where intercepts are in the big theta vector
  rhss <- attr(formula1, "rhs")
  leneqn <- numeric(neqn)
  for(q in 1:neqn){
    leneqn[q] <- length(attr(terms(formula(formula1, lhs = 0, rhs = q)), "term.labels"))
  }
  # in eqn where there are regressors,
  # terms does not count intercept, so we need to add it
  leneqn <- ifelse(leneqn == 1, leneqn, leneqn + 1)
  # print(leneqn)
  # intercept indicies
  # print(cumsum(leneqn)[-c(neqn)]+1)
  int_ind <- c(1, cumsum(leneqn)[-c(neqn)]+1)
  # print(int_ind)
  
  if( is.null(int_ind)) int_ind <- int_inds
  
  
  # ===========================================================
  # search for good initial theta_0
  
  feasible.repeat <- 1000
  myrange <- c(1, 5, 10, 25, 100, 1000)
  
  # ===========================================================
  # step 0: check if can calculate LL at zeros
  
  theta0 <- rep(0, sum(leneqn))
  
  if(!is.null(init)) theta0 <- init
  
  # print(theta0)
  
  # theta0[1] <- NA
  ll0 <- logL(theta0)
  # print(ll0)
  if(is.na(ll0) | ll0 == -Inf){
    cat(paste("\tInitial:           log likelihood = -<Inf> \n", sep = ""), sep = "")
    if(auxilary.reports){
      cat("    searching for feasible values\n", sep = "")
    }
    
    # ===========================================================
    # step 1: finding feasible starting values
    
    myrange <- c(-100, -25, -10, -5, -1, 1, 5, 10, 25, 100, 1000)
    ir <- 1 # ordinal number in "myrange"
    kr <- 1 # ordinal number of the current try of maximum "feasible.repeat"
    ll.feasible <- ll0
    while(is.na(ll.feasible) | ll.feasible == -Inf & kr <= feasible.repeat){
      # cons <- runif(1) * 2 * myrange[ir] + myrange[ir]
      # beta0 <- solve(t(x) %*% x) %*% t(x) %*% rep(cons, n)
      for( qq in int_ind){
        theta0[qq] <- runif(1) * 2 * myrange[ir] + myrange[ir]
      }
      # print(theta0)
      ll.feasible <- logL(theta0)
      ir <- ir + 1
      if(ir > length(myrange)){
        ir <- 1
      }
      kr <- kr + 1
    } # end of while(is.na(ll.feasible) & kr <= feasible.repeat)){
    
    if(is.na(ll.feasible) | ll.feasible == -Inf){
      cat("    could not find feasible values\n", sep = "")
    } else {
      cat(paste("\tFeasible:          log likelihood = ",format(ll.feasible, digits = 13)," \n", sep = ""), sep = "")
      ll0 <- ll.feasible
    }
    
  } else { # end of if(is.na(ll0) | ll0 == -Inf){
    cat(paste("\tInitial:           log likelihood = ",format(ll0, digits = 13)," \n", sep = ""), sep = "")
  }
  
  # logL( theta0 )
  
  
  # ===========================================================
  # step 2: try to improve logL by trying intercepts
  # (exactly 10 tries, not while better is found) randomly
  
  alt.repeat <- length(myrange) * 2 + 1
  
  ir <- 1 # ordinal number in "myrange"
  kr <- 1 # ordinal number of the current try of maximum "alt.repeat"
  ll_alt <- ll0 - 1
  while(kr <= alt.repeat){
    # cons <- runif(1) * 2 * myrange[ir] + myrange[ir]
    # beta1 <- solve(t(x) %*% x) %*% t(x) %*% rep(cons, n)
    theta1 <- theta0
    for( qq in int_ind){
      theta1[qq] <- runif(1) * 2 * myrange[ir] + myrange[ir]
    }
    # print(c(ir,kr,myrange[ir],theta1))
    ll_alt <- logL( theta1 )
    # print(ll_alt)
    # print(c(-1,ll_alt,theta1))
    if(!is.na(ll_alt) & ll_alt != -Inf){
      if(ll_alt > ll0){
        ll0 <- ll_alt
        theta0 <- theta1
      }
    } # end of if(!is.na(ll_alt)){
    ir <- ir + 1
    if(ir > length(myrange)){
      ir <- 1
    }
    kr <- kr + 1
  } # end of while(kr <= alt.repeat){
  
  cat(paste("\tAlternative:       log likelihood = ",format(ll0, digits = 13)," \n", sep = ""), sep = "")
  
  # logL( theta0 )
  
  # ===========================================================
  # step 3: try to improve logL where all coefficients jointly
  # by multiplying them by the same (some) constant c.
  # Since we do not know whether to scale theta up or down,
  # we try both: first down, then up
  
  # print(theta0[1:length(theta0])
  if(auxilary.reports){
    cat("    rescaling the whole vector\n", sep = "")
  }
  
  stheta0 <- 0.5 * theta0
  ll_s <- logL(stheta0)
  # print(c(0,ll_s,theta0[1:length(theta0)]))
  if( !is.na(ll_s) & ll_s != -Inf & ll_s > ll0 ){
    while( !is.na(ll_s) & ll_s != -Inf & ll_s > ll0 ){
      # print(ll_s)
      ll0 <- ll_s
      ll_s <- ll0 - 1
      stheta0 <- 0.5 * stheta0
      ll_s <- logL(stheta0)
      # print(c(1,ll_s,theta0[1:length(theta0)]))
    } # end while
    theta0 <- 2 * stheta0
  } else { # end if( !is.na(ll_s) & ll_s != -Inf & ll_s > ll0 ){ @@@
    stheta0 <- 2 * theta0
    ll_s <- logL(stheta0)
    # print(c(3,ll_s,theta0[1:length(theta0)]))
    if( !is.na(ll_s) & ll_s != -Inf & ll_s > ll0 ){
      while( !is.na(ll_s) & ll_s != -Inf & ll_s > ll0 ){
        ll0 <- ll_s
        ll_s <- ll0 - 1
        stheta0 <- 2 * stheta0
        ll_s <- logL(stheta0)
        # print(c(4,ll_s,theta0[1:length(theta0)]))
      } # end while
      theta0 <- 0.5 * stheta0
      # print(c(5,ll_s,theta0[1:length(theta0)]))
    } else {
      stheta0 <- 0.5 * stheta0
      # print(c(6,ll_s,theta0[1:length(theta0)]))
    }
  } # end else @@@
  cat(paste("\tRescale:           log likelihood = ",format(ll0, digits = 13)," \n", sep = ""), sep = "")
  
  # logL( theta0 )
  
  # ===========================================================
  # step 4: try to improve logL by varying intercepts
  # equation by equation:
  
  # theta11 <- theta0
  if(auxilary.reports){
    cat("    rescaling the vector of slopes: not done\n", sep = "")
    cat("    rescaling the intercepts, eqn by eqn \n", sep = "")
  }
  
  # theta0 <- theta11
  # ww <- seq(M)[-1]	# the column of the theta
  
  ll00 <- ll0
  qq <- 1
  
  repeat{
    
    for( ww in int_ind){
      # print(ww)
      
      stheta0 <- theta0
      stheta0[ww] <- 0.5 * stheta0[ww]
      # ll0 <- logL(theta0)
      ll_s <- logL(stheta0)
      # print(ll_s)
      # print(c(ww-1,1,ll_s,stheta0[1:length(theta0)]))
      if( !is.na(ll_s) & ll_s != -Inf & ll_s > ll0 & abs(ll_s - ll0) > 10^-12){
        # if better, improve it even more
        while( !is.na(ll_s) & ll_s != -Inf & ll_s > ll0 & abs(ll_s - ll0) > 10^-12){
          ll0 <- ll_s
          ll_s <- ll0 - 1
          stheta0[ww] <- 0.5 * stheta0[ww]
          ll_s <- logL(stheta0)
          
          if(report.rescaling){
            cat("    cons_",ww," = ",stheta0[ww],", Case 1a, ll0 = ",format(ll0, digits = 13),", lls = ",format(ll_s, digits = 13),"\n", sep = "")
          }
          # print(c(stheta0[1:length(theta0)]))
        } # end of while
        stheta0[ww] <- 2 * stheta0[ww]
        theta0 <- stheta0
        if( abs( mean( stheta0[ww] ) ) < 10^-8){
          # need to reverse the sign if beta0 gets very small
          stheta0[ww] <- -4 * stheta0[ww]
          ll_s <- logL(stheta0)
          if( !is.na(ll_s) & ll_s != -Inf & ll_s > ll0 & abs(ll_s - ll0) > 10^-12){
            while( !is.na(ll_s) & ll_s != -Inf & ll_s > ll0 & abs(ll_s - ll0) > 10^-12){
              ll0 <- ll_s
              ll_s <- ll0 - 1
              stheta0[ww] <- 2 * stheta0[ww]
              ll_s <- logL(stheta0)
              if(report.rescaling){
                cat("    cons_",ww," = ",stheta0[ww],", Case 1b, ll0 = ",format(ll0, digits = 13),", lls = ",format(ll_s, digits = 13),"\n", sep = "")
              }
              # print(c(stheta0[1:length(theta0)]))
            } # end of while
            stheta0[ww] <- 0.5 * stheta0[ww]
            theta0 <- stheta0
          } else {
            stheta0[ww] <- -0.25 * stheta0[ww]
            theta0 <- stheta0
          } # end else: returned to pre-reserved sign state
        } # end of if( abs( mean( sb0 ) ) < 10^-8){
      } else { # end of if halving was good
        stheta0 <- theta0
        stheta0[ww] <- 2 * stheta0[ww]
        ll_s <- logL(stheta0)
        # print(ll_s)
        # print(c(ww-1,2,ll_s,stheta0))
        # print(!is.na(ll_s) & ll_s != -Inf & ll_s > ll0 & abs(ll_s - ll0) > 10^-12 )
        if( !is.na(ll_s) & ll_s != -Inf & ll_s > ll0 & abs(ll_s - ll0) > 10^-12){
          while( !is.na(ll_s) & ll_s != -Inf & ll_s > ll0 & abs(ll_s - ll0) > 10^-12){
            ll0 <- ll_s
            ll_s <- ll0 - 1
            stheta0[ww] <- 2 * stheta0[ww]
            ll_s <- logL(stheta0)
            if(report.rescaling){
              cat("    cons_",ww," = ",stheta0[ww],", Case 2, ll0 = ",format(ll0, digits = 13),", lls = ",format(ll_s, digits = 13),"\n", sep = "")
            }
            # print(theta0[1:length(theta0])
          } # end of while
          stheta0[ww] <- 0.5 * stheta0[ww]
          theta0 <- stheta0
        } else {# end of if
          stheta0[ww] <- 0.5 * stheta0[ww]
          theta0 <- stheta0
        }
      } # end of else: doubling
    }
    
    if( abs(ll00 - ll0) < 10^-9 | qq > 1){
      break
    }
    
    qq <- qq + 1
    
  }
  
  cat(paste("\tRescale eq:        log likelihood = ",format(ll0, digits = 13)," \n", sep = ""), sep = "")
  
  # logL( theta0 )
  
  return(list(theta0 = theta0, ll0 = ll0))
}

.mlmaximize <- function(theta0, ll, gr = NULL, hess = NULL, alternate = NULL, BHHH = F, level = 0.99, step.back = .Machine$double.eps, reltol = .Machine$double.eps, lmtol = sqrt(.Machine$double.eps), steptol = sqrt(.Machine$double.eps), digits = 4, when.backedup = sqrt(.Machine$double.eps), max.backedup = 11, print.level = 6, only.maximize = FALSE, maxit = 150, sample.size = 100, ...){
  
  theta00 <- theta0
  
  k4 <- length(theta0)
  
  if(print.level >= 6){
    cat("\n=================")
    cat(" Initial values:\n\n", sep = "")
    print(theta0)
  }
  
  #  if(print.level >= 2){
  #   cat("\n=================")
  #   cat(" Maximization:\n\n", sep = "")
  #  }
  
  # step.back = 2^-217
  # print(ll)
  
  ll0 <- ll(theta0, ...)
  # print(ll0)
  ltol <- reltol * (abs(ll0) + reltol)
  typf <- ll0
  theta1 <- theta0
  
  iter <- iter.total <- backedup <- backedups <- wasconcave <- wasconcaves <- 0
  
  if( is.na(ll0) | ll0 == -Inf ){
    if(print.level >= 2){
      cat("Could not compute ll at starting values: trying something else\n")
    }
    iter1 <- backedups
    repeat{
      iter1 <- iter1 + 1
      # theta0 <- theta00 * runif(length(theta0), 0.999, 1.001) # not sure what to do here
      theta0 <- theta00 * rep(0, length(theta0)) # not sure what to do here
      ll0 <- ll( theta0, ... )
      if(!is.na(ll0)) break
      if(iter1 == 55){
        print(iter1)
        # print(ll0)
        stop("it's not gonna happen...")
      }
    }
    # backedups <- iter1
  }
  
  delta1 <- gHg <- s1 <- 1
  h1 <- tryCatch( 2, error = function(e) e )
  cant.invert.hess <- FALSE
  
  if(print.level >= 2){
    cat(paste("Iteration ",formatC(iter, width = 3)," (at starting values):       log likelihood = ",format(ll0, digits = 13),"\n\n", sep = ""), sep = "")
  }
  
  convergence <- -999
  
  repeat{
    iter.total <- iter.total + 1
    # cat("backedup = ",backedup," backedups = ",backedups,"\n", sep = "")
    if(print.level >= 6){
      print(theta0)
    }
    
    # cumulate how many times did it backed-up in a row
    # if(s1 <= when.backedup){    
    if(s1 < when.backedup){
      # backedup <- backedup + 1
    } else {
      # backedup <- backedup
      backedup <- 0
    }
    # cat.print(s1)
    # cat.print(when.backedup)
    # cat.print(backedup)
    
    # print(s1)
    # cumulate how many times was concave
    if( inherits(h1, "error") ){
      wasconcave <- wasconcave + 1
    } else {
      wasconcave <- 0
    }
    
    # try different values if was concave more than @@@ times
    if(wasconcave == max.backedup){
      # start over
      if(print.level >= 2){
        cat("Not concave ",max.backedup," times in a row: trying something else (not concave ",wasconcaves+1," times in total)\n")
      }
      iter <- wasconcave <- backedup <- 0
      wasconcaves <- wasconcaves + 1
      theta0 <- theta0*runif(length(theta0), 0.999, 1.001) # not sure what to do here
      ll0 <- ll( theta0, ... )
      # if(print.theta) print(theta0)
      if(print.level >= 2){
        cat(paste("Iteration ",formatC(iter, width = 3),"  (at slightly perturbed starting values):       log likelihood = ",format(ll0, digits = 13),"\n\n", sep = ""), sep = "")
      }
    }
    
    # try different values if backed-up more than @@@ times
    if(backedup == max.backedup){
      # start over
      if(print.level >= 2){
        cat("Backed-up ",max.backedup," times in a row: trying something else (backup-up ",backedups+1," times in total)\n", sep = "")
      }
      iter <- backedup <- wasconcave <- 0
      backedups <- backedups + 1
      wasconcaves <- wasconcaves + 1
      theta0 <- theta0*runif(length(theta0), 0.999, 1.001) # not sure what to do here
      ll0 <- ll( theta0, ... )
      if(print.level >= 6){
        print(theta0)
      }
      if(print.level >= 2){
        cat(paste("Iteration ",formatC(iter, width = 3)," (at slightly perturbed starting values):       log likelihood = ",format(ll0, digits = 13),"\n\n", sep = ""), sep = "")
      }
    }
    
    # see if it calculated ll
    if( is.na(ll0) | ll0 == -Inf | ll0 == Inf | ll0 == 0 ){
      if(print.level >= 2){
        cat("Could not compute ll: trying something else\n")
      }
      iter1 <- backedups
      repeat{
        iter1 <- iter1 + 1
        # theta0 <- c( cons0, beta0, mu = 0, eta = 0, lnsv2 = -1*iter1/2, lnsu2 = -1*iter1/2)
        theta0 <- theta00*runif(length(theta0), 0.999, 1.001) # not sure what to do here
        ll0 <- ll( theta0, ... )
        if(!is.na(ll0) & ll0 != 0) break
        if(iter1 == 15){
          stop("it's not gonna happen... could not compute at initial and find feasible values")
        }
      }
      iter <- backedup <- 0
      backedups <- iter1
      if(print.level >= 2){
        cat(paste("Iteration ",formatC(iter, width = 3)," (at slightly perturbed starting values):       log likelihood = ",format(ll0, digits = 13),"\n\n", sep = ""), sep = "")
      }
    }
    if(backedups == 500){
      stop("it's not gonna happen... backuped 5 times")
    }
    if(wasconcaves == 500){
      stop("it's not gonna happen... not concave 5 times")
    }
    iter <- iter + 1
    delta3 <- 1
    # step 2: direction vector
    
    # previous theta
    
    if(iter.total > 1) theta1 <- theta0 - s1 * d0
    
    # BHHH (faster, but different):
    # The Hessian is approximated as the negative
    # of the sum of the outer products of the gradients
    # of individual observations,
    # -t(gradient) %*% gradient = - crossprod( gradient )
    g1 <- gr(theta0, ...)
    # print(g1)
    if(!is.null(alternate)) BHHH <- floor(iter/alternate) %% 2 == 1
    h0 <- hess(theta0,  ...)
    # print(h0)
    # check if negative definite
    eigen1 <- eigen( h0 )
    # eigen.tol <- k4 * max(abs(eigen1$values)) * .Machine$double.eps # this is for positive definiteness
    eigen.val <- eigen1$values#ifelse(eigen1$values < .Machine$double.eps^.1, -1e-4, eigen1$values)
    # hess.pos.def <- sum(eigen1$values > eigen.tol) == k4
    hess.neg.def <- !any(eigen.val >= 0)
    # 1. replace negative with small ones
    # eigen.val <- ifelse(eigen1$values < 0, .0001, eigen.val)
    # 2. replace negative with absolut values
    eigen.val <- abs(eigen1$values)
    # print(hess.neg.def)
    # make it negative definite if it is not already
    if(!hess.neg.def){
      h0_ <- matrix(0, k4, k4)
      # eigen1 <- eigen( h0 )
      for( i in seq_len( k4 ) ){
        # h0_ <- h0_ - abs(eigen1$values[i]) * eigen1$vectors[,i] %*% t(eigen1$vectors[,i])
        h0_ <- h0_ - eigen.val[i] * eigen1$vectors[,i] %*% t(eigen1$vectors[,i])
      }
      h0.old <- h0
      h0 <- h0_
    }
    # print( is.negative.definite(h0) )
    # remember hessian and negative of its inverse from previous iter that could have been inverted
    if( !cant.invert.hess ){
      h0_previous <- h0
      h1_previous <- h1
    }
    # h0 <- h0.old
    # easier to invert positive definite matrix
    h1 <- tryCatch( qr.solve(-h0, tol = 1e-10), error = function(e) e )
    # check if it can be inverted
    cant.invert.hess <- FALSE
    cant.invert.hess <- inherits(h1, "error")
    if( cant.invert.hess ){
      # print(h1)
      if(print.level >= 2){
        # cat(paste("cannot invert Hessian, using eigenvalues\n", sep = ""), sep = "")
      }
      # this was just to get the uninvertable hessian
      # return(list(hess = h0, grad = g1))
      # stop("this")
      # @14@ this
      eig1 <- eigen( -h0_previous )
      d0 <- rep(0, length(theta0))
      # eig2 <- ifelse(eig1$values < eps1, 1, eig1$values)
      for (i in 1:length(eig1$values)){
        d0 <- d0 + (g1%*%eig1$vectors[,i])*eig1$vectors[,i] / eig1$values[i]
      }
      # @14@ could be done easier
      # d0 <- qr.solve(-h0, g1, tol = 1e-10)
      gHg <- sum( g1 * d0)
      # in the part of the ortogonal subspace where the eigenvalues
      # are negative or small positive numbers, use steepest ascent
      # in other subspace use NR step
      # d0 <- ifelse(eigen(-h0, only.values = TRUE)$values < reltol, g1, d0)
      gg <- sqrt( crossprod(g1) )
      gHg <- gg
      # d0 <- g1
      # d0
    } else {
      d0 <- as.vector( h1 %*% g1 )
      gg <- sqrt( crossprod(g1) )
      # h1.old <- solve(-h0.old)
      gHg <- as.vector( t(g1) %*% h1 %*% g1 )
    }
    # gg_scaled <- gg * max( crossprod(theta0), crossprod(theta1) ) / max( abs(ll0), abs(typf))
    # theta_rel <- max( abs(theta0 - theta1) / apply( cbind( abs(theta0),abs(theta1) ), 1, max) )
    theta_rel <- max( abs(theta0 - theta1) / (abs(theta1)+1) )
    
    
    # begin stopping criteria calculated using new values of g1 and h1
    if(s1 > when.backedup & delta1 != 17.17){ # if(s1 > when.backedup*10^-100 & !cant.invert.hess){
      if(abs(gHg) < lmtol & iter.total > 1){
        convergence <- 0
        if(print.level >= 2){
          cat("\nConvergence given g inv(H) g' = ",abs(gHg)," < lmtol\n", sep = "")
        }
        break
      }
      if(theta_rel < steptol & iter.total > 2){
        # print(theta_rel)
        convergence <- 0
        if(print.level >= 2){
          cat("\nConvergence given relative change in parameters = ",theta_rel," < steptol\n", sep = "")
        }
        break
      }
    }
    # end stopping criteria
    # use steepest ascent when backed-up
    if(s1 <= when.backedup*10^-0){
      # eig1 <- eigen( -h0 )
      d0 <- g1*runif(length(g1), min = 0.95, max = 1.05)
      # d0 <- ifelse(eig1$values < reltol, g1, d0)
      # theta0 <- theta0 - 1 * d0
    }
    # print(d0)
    # step 3: new guess
    # a: s = 1
    # b: funct(theta0 + d0) > funct(theta0)
    s1 <- 1
    
    s1 <- sum(theta0^2)/abs(gHg)
    
    theta1 <- theta0 + s1 * d0
    # print(12)
    # print(theta1)
    ll1 <- ll( theta1, ... )
    # print(13)
    delta2 <- ll1 - ll0
    flag <- (!is.na(delta2) & delta2 != -Inf & delta2 > 0)
    # begin Cases
    if( flag ){
      # begin Case 1: f(theta1) > f(theta0)
      ll.temp <- ll0
      # check if s1 = 2, 3, ... increases f even more
      while( flag ){
        if(print.level >= 6){
          cat(paste("\t\tCase 1: s = ",s1,", delta = ",delta2,"\n", sep = ""), sep = "")
        }
        ll0 <- ll1
        s1 <- s1 + 1
        theta1 <- theta0 + s1 * d0
        ll1 <- ll( theta1, ... )
        delta2 <- ll1 - ll0
        flag <- (!is.na(delta2) & delta2 != -Inf & delta2 > 0)
      }
      # overall delta
      delta1 <- ll0 - ll.temp
      delta_rel <- abs(delta1 / ll.temp)
      # print(delta_rel)
      s1 <- s1 - 1
      # overwrite the values
      theta0 <- theta0 + s1 * d0
      # end Case 1: f(theta1) > f(theta0)
    } else {
      # begin Case 2: f(theta1) < f(theta0)
      # check only if s1=1/2 increases f
      s1 <- 0.5
      
      s1 <- 0.5 * sum(theta0^2)/abs(gHg)
      
      theta1 <- theta0 + s1 * d0
      ll1 <- ll( theta1, ... )
      # cat(" ll1 = ",ll1,", ll0 = ",ll0,", s = ",s1,"\n", sep = "")
      delta2 <- ll1 - ll0
      if(print.level >= 6){
        cat(paste("\t\tCase 2: s = ",s1,", delta = ",delta2,"\n", sep = ""), sep = "")
      }
      flag2 <- (!is.na(delta2) & delta2 != -Inf & delta2 > 0)
      # end Case 2: f(theta1) < f(theta0)
      if( flag2 ){
        # begin Case 2a: f(theta1) > f(theta0)
        ll.temp <- ll0
        # check if s1=1/2^2,1/2^3,... increases f even more
        while( flag2 ){
          if(print.level >= 6){
            cat(paste("\t\t\tCase 2a: s = ",s1,", delta = ",delta2,"\n", sep = ""), sep = "")
          }
          ll0 <- ll1
          s1 <- 0.5 * s1
          theta1 <- theta0 + s1 * d0
          ll1 <- ll( theta1, ... )
          # cat(" ll1 = ",ll1,", ll0 = ",ll0,", s = ",s1,"\n", sep = "")
          delta2 <- ll1 - ll0
          flag2 <- (!is.na(delta2) & delta2 != -Inf & delta2 > ltol)
        }
        # overall delta
        delta1 <- ll0 - ll.temp
        delta_rel <- abs(delta1 / ll.temp)
        # print(delta_rel)
        s1 <- 2 * s1
        # overwrite the values
        theta0 <- theta0 + s1 * d0
        # end Case 2a: f(theta1) > f(theta0)
      } else {
        # begin Case 2b: f(theta1) < f(theta0)
        ll.temp <- ll0
        # try s1=1/2^2,1/2^3,... so that f(theta1) > f(theta0)
        while ( !flag2 & s1 > step.back ){
          s1 <- 0.5 * s1
          theta1 <- theta0 + s1 * d0
          ll1 <- ll( theta1, ... )
          delta2 <- ll1 - ll0
          if(print.level >= 6){
            cat(paste("\t\t\tCase 2b: s = ",s1,", delta = ",delta2,"\n", sep = ""), sep = "")
          }
          flag2 <- (!is.na(delta2) & delta2 != -Inf & delta2 > 0)
        }
        if( !flag2 | s1 < step.back ){
          # stop("provide different starting values")
          delta1 <- 17.17
        } else {
          # overwrite the values
          delta1 <- delta2
          delta_rel <- abs(delta1 / ll.temp)
          ll0 <- ll1
          theta0 <- theta0 + s1 * d0
        }
        # end Case 2b: f(theta1) < f(theta0)
      }
    }
    
    
    if(print.level >= 2){
      if( cant.invert.hess ){
        cat(paste("Iteration ",formatC(iter, width = 3)," (hessian is ",ifelse(BHHH, "BHHH", "provided"),", ",formatC(iter.total, width = 3)," in total):   log likelihood = ",format(ll0, digits = 13)," (not concave)\n", sep = ""), sep = "")
      } else if (s1 <= when.backedup) {
        cat(paste("Iteration ",formatC(iter, width = 3)," (hessian is ",ifelse(BHHH, "BHHH", "provided"),", ",formatC(iter.total, width = 3)," in total):   log likelihood = ",format(ll0, digits = 13)," (backed up)\n", sep = ""), sep = "")
      } else {
        cat(paste("Iteration ",formatC(iter, width = 3)," (hessian is ",ifelse(BHHH, "BHHH", "provided"),", ",formatC(iter.total, width = 3)," in total):   log likelihood = ",format(ll0, digits = 13),"\n", sep = ""), sep = "")
      }
    }
    
    # printing criteria
    if(print.level >= 5){
      if( cant.invert.hess ){
        cat(paste(" (in iter ",formatC(iter, width = 3),": delta = ",format(delta1, digits = 6),"; s = ",format(s1, digits = 6),"; quasi-gHg = ",format(gHg, digits = 6),"; theta_rel_change = ",format(theta_rel, digits = 6),")\n\n", sep = ""), sep = "")
        if(print.level >= 5.5){
          print(theta0)
          cat("\n")
        }
      } else {
        cat(paste(" (in iter ",formatC(iter, width = 3),": delta = ",format(delta1, digits = 6),"; s = ",format(s1, digits = 6),"; gHg = ",format(gHg, digits = 6),"; theta_rel_change = ",format(theta_rel, digits = 6),")\n\n", sep = ""), sep = "")
        if(print.level >= 5.5){
          print(theta0)
          cat("\n")
        }
      }
    }
    # print(s1)
    if(s1 > when.backedup & !cant.invert.hess){ # if(s1 > when.backedup^2 & !cant.invert.hess){
      # ltol <- reltol * (abs(ll0) + reltol)
      # print(cant.invert.hess)
      if(delta1 > 0 & !is.na(delta_rel) & delta_rel < ltol*1e-5 & iter.total > 1){
        convergence <- 0
        if(print.level >= 2){
          cat("\nConvergence given relative change in log likelihood = ",delta_rel," < ltol\n", sep = "")
        }
        break
      }
    }
    if(iter.total > maxit){
      convergence <- 1
      cat("\n Maximum number of iterations (",maxit,") reached without convergence\n", sep = "")
      cat("convergence <- 1\n")
      break
    }
  } # end repeat
  
  if( !only.maximize & !cant.invert.hess){
    names(ll0) <- NULL
    colnames(h1) <- rownames(h1) <- names(g1) <- names(theta0)
    
    # sqrt(crossprod(g1))
    
    b0 <- theta0
    sd0 <- sqrt( diag( h1 ) )
    t0 <- b0 / sd0
    p0 <- pt(abs(t0), sample.size-length(b0), lower.tail = FALSE) * 2
    t10 <- qt((1-0.01*level)/2, sample.size-length(b0), lower.tail = FALSE)
    t17 <- cbind( b0, sd0, t0, p0, b0 - t10*sd0, b0 + t10*sd0)
    # t17 <- cbind( b0, sd0, t0, p0)
    colnames(t17) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)", paste("",level,"_CI_LB", sep = ""), paste("",level,"_CI_UB", sep = ""))
    # colnames(t17) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    # t17
    
    if(print.level >= 2){
      cat(paste("\nFinal log likelihood = ",format(ll0, digits = 13),"\n\n", sep = ""), sep = "")
    }
    # cat(paste("Stoc. frontier normal/",distribution,"\n", sep = ""), sep = "")
    if(print.level >= 5.5){
      cat("\nCoefficients:\n\n", sep = "")
      printCoefmat(t17[,1:4], digits = digits)
    }
    
    return(list(par = theta0, table = t17, gradient = g1, vcov = h1, ll = ll0, gg = gg, gHg = gHg, theta_rel_ch = theta_rel, convergence = convergence))
  } else {
    return(list(par = theta0, gradient = g1, vcov = h1, ll = ll0, gg = gg, gHg = gHg, theta_rel_ch = theta_rel, convergence = convergence))
  }
  
  
}

gradient1 <- function(func, at, ...){
  K <- length(at)
  hs <- (.Machine$double.eps)^(1/3) * ifelse(abs(at) > 1, abs(at), 1)
  grad2 <- numeric(K)
  for(i in seq(K)){
    # cat("   Current derivative is ",i," (",i," of ",K,")\n", sep = "")
    zeros.i <- numeric(K); zeros.i[i] <- hs[i] / 2
    grad2[i] <- (func( at + zeros.i, ... ) - func( at - zeros.i, ... ) ) / hs[i]
  }
  return(grad2)
  # cat("   \n", sep = "")
}

gradient1.by.i <- function(func, at, ...){
  K <- length(at)
  hs <- (.Machine$double.eps)^(1/3) * ifelse(abs(at) > 1, abs(at), 1)
  grad2 <- NULL
  for(i in seq(K)){
    # cat("   Current derivative is ",i," (",i," of ",K,")\n", sep = "")
    zeros.i <- numeric(K); zeros.i[i] <- hs[i] / 2
    tymch <- (func( at + zeros.i, ... ) - func( at - zeros.i, ... ) ) / hs[i]
    grad2 <- cbind(grad2, tymch)
    # if(i < 2){
    #   cat.print(tymch)
    # }
  }
  colnames(grad2) <- names(at)
  # cat.print(grad2[1:10,])
  return(grad2)
  # cat("   \n", sep = "")
}

gradient1.mp <- function(func, at, ...){
  # dots <- list(...)
  # cat.print(names( list(...) ) )
  K <- length(at)
  hs <- (.Machine$double.eps)^(1/3) * ifelse(abs(at) > 1, abs(at), 1)
  # cat.print(hs)
  # cat.print( list(...)$constr.ort)
  constr    <- !is.null( list(...)$constr.ort) & sum( list(...)$constr.ort) > 0
  if( constr ){
    constr.in <- which( list(...)$constr.ort)
    hs[ list(...)$constr.ort] <- 0
  }
  # cat.print(hs)
  ats <- list()
  for (k in 1:K) {
    # cat.print(k)
    # ats[[ (k-1)*2 + 1]]    <- ats[[ k*2        ]]    <- at
    # ats[[ (k-1)*2 + 1]][k] <- ats[[ (k-1)*2 + 1]][k] + hs[k] / 2
    # ats[[ k*2        ]][k] <- ats[[ k*2        ]][k] - hs[k] / 2    
    ats[[ 2*k - 1 ]]    <- ats[[ 2*k     ]]                <- at
    ats[[ 2*k - 1 ]][k] <- ats[[ 2*k - 1 ]][k] + hs[k] / 2
    ats[[ 2*k     ]][k] <- ats[[ 2*k     ]][k] - hs[k] / 2
  }
  # cat.print(ats)
  tymch1 <- foreach(i = seq_len(2*K), .combine = c) %dopar%{
    func( ats[[i]], ... )
  }
  # cat.print(tymch1)
  grad2 <- numeric(K)
  for (k in 1:K) {
    grad2[k] <- ( tymch1[[ 2*k - 1 ]] - tymch1[[2*k]] ) / hs[k]
  }
  if( constr ){
    grad2[constr.in] <- 0
  }
  return(grad2)
  # cat("   \n", sep = "")
}

hessian1 <- function(func, at, print = FALSE, ...){
  # cat.print(at)
  K <- length(at)
  hs <- (.Machine$double.eps)^(1/4) * ifelse(abs(at) > 1, abs(at), 1)
  hess2 <- matrix(0, nrow = K, ncol = K)
  fz <- func( at, ... )
  for(i in seq(K)){
    zeros.i <- numeric(K); zeros.i[i] <- hs[i]
    hess2[i,i] <- (	func( at + zeros.i, ... ) + func( at - zeros.i, ... ) - 2*fz ) / hs[i] / hs[i]
    # print(hess2[i,i])
    # if(print) cat("before loop over j", i,"\n")
    zeros.i[i] <- hs[i] / 2
    for(j in seq(K)){
      if(i <= j){
        if(print) cat("   i = ",i,"; j = ",j," (out of ",K,"x",K,")\n", sep = "")
        # if(print) cat("in the loop over j", i, j,"\n")
        zeros.j <- numeric(K); zeros.j[j] <- hs[j] / 2
        hess2[i,j] <- (	func( at + zeros.j + zeros.i, ... ) + func( at - zeros.j - zeros.i, ... ) - func( at + zeros.j - zeros.i, ... ) - func( at - zeros.j + zeros.i, ... ) ) / hs[i] / hs[j]
        # print(hess2[i,j])
      }
    }
  }
  # cat("   \n", sep = "")
  hess2 <- hess2 + t(hess2) - diag(diag(hess2))
  return(hess2)
}

hessian1.mp <- function(func, at, ...){
  
  K         <- length(at)
  hs        <- (.Machine$double.eps)^(1/4) * ifelse(abs(at) > 1, abs(at), 1)
  # cat.print(hs)
  # cat.print(constr.ort)
  constr    <- !is.null( list(...)$constr.ort) & sum( list(...)$constr.ort) > 0
  if( constr ){
    constr.in <- which( list(...)$constr.ort)
    hs[ list(...)$constr.ort] <- 0
  }
  # cat.print(hs)
  # we have K*(K+1)/2 elements to fill
  # for each of these elements we need 4 elements
  # our at is 2*K*(K+1)
  # p ~ plus
  # m ~ minus
  ats <- list()
  for(i in 1:K){
    for(j in 1:K){
      if(i <= j){
        myind <- 4*( (i-1)*K - (i-1)*(i-2)/2 ) + 4*j - (i-1)*4
        # if(i >= 2) myind <- myind - 4
        # cat("i = ",i,", j = ",j,", myind = ",myind, " element 1",myind-3,"\n")
        # at.p.j.p.i
        ats[[ myind - 3 ]]    <- at
        ats[[ myind - 3 ]][j] <- ats[[ myind - 3 ]][j] + hs[j] / 2
        ats[[ myind - 3 ]][i] <- ats[[ myind - 3 ]][i] + hs[i] / 2
        # at.m.j.m.i
        ats[[ myind - 2 ]]    <- at
        ats[[ myind - 2 ]][j] <- ats[[ myind - 2 ]][j] - hs[j] / 2
        ats[[ myind - 2 ]][i] <- ats[[ myind - 2 ]][i] - hs[i] / 2
        # at.p.j.m.i
        ats[[ myind - 1 ]]    <- at
        ats[[ myind - 1 ]][j] <- ats[[ myind - 1 ]][j] + hs[j] / 2
        ats[[ myind - 1 ]][i] <- ats[[ myind - 1 ]][i] - hs[i] / 2
        # at.m.j.p.i
        ats[[ myind - 0 ]]    <- at
        ats[[ myind - 0 ]][j] <- ats[[ myind - 0 ]][j] - hs[j] / 2
        ats[[ myind - 0 ]][i] <- ats[[ myind - 0 ]][i] + hs[i] / 2
      }
    }
  }
  # cat.print(ats)
  if( constr ){
    tymch1 <- foreach(k = seq_len(length(ats))) %dopar%{
      func( ats[[k]], ... )
    }
  } else {
    tymch1 <- foreach(k = seq_len(length(ats))) %dopar%{
      func( ats[[k]], ... )
    }
  }
  # cat.print(tymch1)
  # collect and fill
  hess2 <- matrix(0, nrow = K, ncol = K)
  for(i in 1:K){
    for(j in 1:K){
      if(i <= j){
        myind <- 4*( (i-1)*K - (i-1)*(i-2)/2 ) + 4*j - (i-1)*4
        hess2[i,j] <- hess2[j,i] <- 
          (tymch1[[ myind - 3 ]] + tymch1[[ myind - 2 ]] -
             tymch1[[ myind - 1 ]] - tymch1[[ myind - 0 ]]) / hs[i] / hs[j]
      }
    }
  }
  if( constr ){
    hess2[constr.in,] <- hess2[,constr.in] <- 0
  }
  # print(hess2)
  # hess2 <- hess2 + t(hess2) - diag(diag(hess2))
  return(hess2)
}

hessian2 <- function(funcg, at, ...){
  K <- length(at)
  hs <- (.Machine$double.eps)^(1/3) * ifelse(abs(at) > 1, abs(at), 1)
  hess0 <- matrix(0, nrow = K, ncol = K)
  colnames(hess0) <- rownames(hess0) <- names(at)
  for(i in seq(K)){
    zeros.i <- numeric(K)
    zeros.i[i] <- hs[i] / 2
    hess0[,i] <- (funcg( at + zeros.i, ... ) - funcg( at - zeros.i, ... ) ) / hs[i]
  }
  hess2 <- ( hess0 + t(hess0) ) / 2
  return(hess2)
}





#' Pretty-print coefficient matrix (internal helper)
#'
#' @keywords internal
#' @noRd
coef.mat.print <- function(coeffs, digits = 4, large = 999) {
  if (!is.matrix(coeffs)) coeffs <- as.matrix(coeffs)
  
  big_or_small <- abs(coeffs[, , drop = FALSE]) > large |
    abs(coeffs[, , drop = FALSE]) < 10^(-digits)
  
  out <- coeffs
  out[big_or_small] <- formatC(coeffs[big_or_small, drop = TRUE],
                               digits = 1, format = "e",
                               width = digits + 5)
  out[!big_or_small] <- formatC(coeffs[!big_or_small, drop = TRUE],
                                digits = digits, format = "f",
                                width = digits + 5)
  out
}

# coef.mat.print <- function(coeffs, digits = 4, large = 999){
#   if(!is.matrix(coeffs)) coeffs <- as.matrix(coeffs)
#   ifelse(abs(coeffs[,, drop = FALSE]) > large | abs(coeffs[,, drop = FALSE]) < 10^-(digits), formatC(coeffs[,, drop = FALSE], digits = 1, format = "e", width = digits+5), formatC(coeffs[,, drop = FALSE], digits = digits, format = "f", width = digits+5))
# }





#

# Stata: M.pdf page 698
# C_concave: -H is positive semidefinite
.is.positive.semi.definite <- function (x, tol = 1e-08){
  if (!is.numeric(x))
    stop("argument x is not a numeric matrix")
  eigenvalues <- eigen(x, only.values = TRUE)$values
  n <- nrow(x)
  for (i in 1:n) {
    if (abs(eigenvalues[i]) < tol) {
      eigenvalues[i] <- 0
    }
  }
  if (any(eigenvalues < 0)) {
    return(FALSE)
  }
  return(TRUE)
}

.is.negative.definite <- function (x, tol = 1e-16)
{
  # if (!is.square.matrix(x))
  # stop("argument x is not a square matrix")
  # if (!is.symmetric.matrix(x))
  # stop("argument x is not a symmetric matrix")
  # if (!is.numeric(x))
  # stop("argument x is not a numeric matrix")
  eigenvalues <- eigen(x, only.values = TRUE)$values
  # cat("Eigenvalues\n")
  # print(eigenvalues)
  n <- nrow(x)
  for (i in 1:n) {
    if (abs(eigenvalues[i]) < tol) {
      eigenvalues[i] <- 0
    }
  }
  if (any(eigenvalues >= 0)) {
    return(FALSE)
  }
  return(TRUE)
}

.make.neg.def <- function(hess){
  # K4 <- ncol(hess)
  # h0_ <- matrix(0, K4, K4)
  eigen1 <- eigen( hess )
  # eigen.val <- abs(eigen1$values)
  #   for( i in seq_len( K4 ) ){
  #     h0_ <- h0_ - abs(eigen1$values[i]) * eigen1$vectors[,i] %*% t(eigen1$vectors[,i])
  #     # h0_ <- h0_ - eigen.val[i] * eigen1$vectors[,i] %*% t(eigen1$vectors[,i])
  #   }
  return(-crossprod(t(eigen1$vectors)*abs(eigen1$values), t(eigen1$vectors)))
}

.make.invertible <- function(hess){
  K <- ncol(hess)
  ok <- FALSE
  adj <- sqrt(.Machine$double.eps)
  # tr0 <- sum( 1/eigen(H,only.values = T)$values )
  i <- 0
  while(!ok & i < 1000){
    i <- 1 + i
    hess <- hess + adj*diag(rep(1, K))
    mytry <- tryCatch( solve(-hess), error = function(e) e )
    ok <- !inherits( mytry, "error")
    # print(mytry)
    # ok <- all(diag(H) != 0)
    adj <- adj * 2
    # print(adj)
  }
  # print(c(i, adj))
  return(hess)
}

.make.zero.one <- function(qq){
  for (i in 1:length(qq)) {
    qqq <- qq[i]
    if(abs(qqq) > 1){
      while(abs(qqq) > 1){
        qqq <- qqq / 10
      }
      qq[i] <- qqq
    }
  }
  return(qq)
}

# http://statsadventure.blogspot.de/2011/12/non-pd-matrices-in-r-cont.html
.nearPSD <- function(c){
  n = dim(c)[1]
  e = eigen(c,sym=TRUE)
  val = e$values * (e$values > 0)
  vec = e$vectors
  T = sqrt(diag( as.vector(1/(vec^2 %*% val)),n,n))
  B = T %*% vec %*% diag(as.vector(sqrt(val)),n,n)
  out = B %*% t(B)
  return(out)
}

.is.negative.definite <- function (x, tol = 1e-16)
{
  # if (!is.square.matrix(x))
  # stop("argument x is not a square matrix")
  # if (!is.symmetric.matrix(x))
  # stop("argument x is not a symmetric matrix")
  # if (!is.numeric(x))
  # stop("argument x is not a numeric matrix")
  eigenvalues <- eigen(x, only.values = TRUE)$values
  # cat("Eigenvalues\n")
  # print(eigenvalues)
  n <- nrow(x)
  for (i in 1:n) {
    if (abs(eigenvalues[i]) < tol) {
      eigenvalues[i] <- 0
    }
  }
  if (any(eigenvalues >= 0)) {
    return(FALSE)
  }
  return(TRUE)
}

.make.neg.def <- function(hess){
  # K4 <- ncol(hess)
  # h0_ <- matrix(0, K4, K4)
  eigen1 <- eigen( hess )
  # eigen.val <- abs(eigen1$values)
  #   for( i in seq_len( K4 ) ){
  #     h0_ <- h0_ - abs(eigen1$values[i]) * eigen1$vectors[,i] %*% t(eigen1$vectors[,i])
  #     # h0_ <- h0_ - eigen.val[i] * eigen1$vectors[,i] %*% t(eigen1$vectors[,i])
  #   }
  return(-crossprod(t(eigen1$vectors)*abs(eigen1$values), t(eigen1$vectors)))
}

.make.invertible <- function(hess){
  K <- ncol(hess)
  ok <- FALSE
  adj <- sqrt(.Machine$double.eps)
  # tr0 <- sum( 1/eigen(H,only.values = T)$values )
  i <- 0
  while(!ok & i < 1000){
    i <- 1 + i
    hess <- hess + adj*diag(rep(1, K))
    mytry <- tryCatch( solve(-hess), error = function(e) e )
    ok <- !inherits( mytry, "error")
    # print(mytry)
    # ok <- all(diag(H) != 0)
    adj <- adj * 2
    # print(adj)
  }
  # print(c(i, adj))
  return(hess)
}

.make.zero.one <- function(qq){
  for (i in 1:length(qq)) {
    qqq <- qq[i]
    if(abs(qqq) > 1){
      while(abs(qqq) > 1){
        qqq <- qqq / 10
      }
      qq[i] <- qqq
    }
  }
  return(qq)
}

cat.print <- function(x, name = NULL){
  if(is.null(name) | !is.character(name)){
    cat(" ",deparse(substitute(x)),":\n", sep = "")
  } else {
    cat(" ",name,":\n", sep = "")
  }
  
  print(x)
}