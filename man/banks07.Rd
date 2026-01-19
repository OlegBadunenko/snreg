\name{banks07}
\docType{data}
\alias{banks07}

\title{U.S. Commercial Banks Data}

\description{
\code{banks07} data frame contains selected variables for 500 (randomly sampled from around 5000) U.S. commercial banks from data of Koetter et al. (2012) for year 2007.  This data are used for illustrution purposes and are not suitable for research purposes. 
}

\usage{data(banks06_07)}

\format{
 This data frame contains the following variables (columns):
  \describe{
  \item{\code{year}}{Year.}
  \item{\code{id}}{Entity ID.}
  \item{\code{TA}}{Gross total assets.}
  \item{\code{LLP}}{Loan loss provisions.}
  \item{\code{Y1}}{Total securities (in thousands of US dollars).}
  \item{\code{Y2}}{Total loans and leases (in thousands of US dollars).}
  \item{\code{W1}}{Cost of fixed assets divided by the cost of borrowed funds.}
  \item{\code{W2}}{Cost of labor (in thousands of US dollars) divided by the cost of borrowed funds.}
  \item{\code{ER}}{Gross total equity to gross total assets ratio.}
  \item{\code{TC}}{Total operating cost.}
  \item{\code{LA}}{Total loans and leases to gross total assets ratio.}}
}

\details{
 The data were sampled and generated as shown in section "Examples".
 
 }
 
\examples{

\dontrun{

data(banks07, package = 'snreg')
dat1 <- banks07
head(dat1)
dat1$mysample <- TRUE

myprod <- FALSE

spe.tl <- 
  log(TC/W3) ~ (log(Y1) + log(Y2) + log(W1/W3) + log(W2/W3))^2 +
  I(0.5*log(Y1)^2) + I(0.5*log(Y2)^2) + 
  I(0.5*log(W1/W3)^2) + I(0.5*log(W2/W3)^2)

# N0 ---------------------------------------

m.N0 <- 
  snreg::lm.mle  (
    formula = spe.tl, data = dat1, subset = mysample,
    ln.var.v  = NULL,
    technique = c('bfgs')
  )

m.N0$value
m.N0$coef
sqrt(diag(m.N0$vcov))
cbind(coef = m.N0$coef, z = m.N0$coef/sqrt(diag(m.N0$vcov)))


# SN0 --------------------------------------

formSK <- NULL
formSV <- NULL


init.skY <- seq(0,3,.07)
init.skY[1] <- 0.01
init.skY <- sort(c(-init.skY,init.skY))

tymch1 <- NULL

for (in.sk in init.skY) {
  m0.sn.tl <- 
    snreg::snreg  (
      formula = spe.tl, data = dat1, subset = mysample,
      start.sk  = in.sk,
      ln.var.v  = formSV,
      skew.v    = formSK,
      start.val = NULL,
      technique = c('bfgs'),
      lmtol     = 1e-5,  
      reltol    = 1e-15, 
      maxit     = 399, 
      report    = 1,
      trace     = 1,
      print.level = 4, 
      digits    = 4)
  tymch1 <- rbind(tymch1, c(in.sk, m0.sn.tl$value))
}

tymch1[order(tymch1[,2]),]

m.SN0 <- 
  snreg::snreg  (
    formula = spe.tl, data = dat1, subset = mysample,
    start.sk  = tymch1[order(tymch1[,2]),][ nrow(tymch1),1],
    ln.var.v  = formSV,
    skew.v    = formSK,
    start.val = NULL,
    technique = c('nr')
  )

m.SN0$ll
m.SN0$coef
sqrt(diag(m.SN0$vcov))
cbind(coef = m.SN0$coef, z = m.SN0$coef/sqrt(diag(m.SN0$vcov)))


# SF1 --------------------------------------


formSK <- NULL
formSU <- NULL
formSV <- NULL

m.SF1 <- 
  snreg::snsf  (
    formula = spe.tl, 
    data = dat1, subset = mysample,
    distribution = "e",
    prod      = myprod,
    start.sk  = -1,
    ln.var.u  = formSU,
    ln.var.v  = formSV,
    skew.v    = formSK,
    mean.u    = NULL,
    start.val = NULL,
    technique = c('bfgs')
  )

m.SF1$ll
m.SF1$coef
sqrt(diag(m.SF1$vcov))
cbind(coef = m.SF1$coef, z = m.SF1$coef/sqrt(diag(m.SF1$vcov)))




# SF2 --------------------------------------


formSK <- ~ SDROA
formSV <- ~ log(TA)
formSU <- ~ scope

m.SF2 <- 
  snreg::snsf  (
    formula = spe.tl, data = dat1, subset = mysample == TRUE,
    distribution = "e",
    prod      = myprod,
    start.sk  = -0.5,
    ln.var.u  = formSU,
    ln.var.v  = formSV,
    skew.v    = formSK,
    mean.u    = NULL,
    start.val = NULL,
    technique = c('bfgs'),
    # vcetype   = 'opg',
    lmtol     = 1e-5,  
    reltol    = 1e-15, 
    maxit     = 999, 
    report    = 1,
    trace     = 1,
    print.level = 5, 
    digits    = 4,
    threads   = 1)


m.SF2$ll
m.SF2$coef
sqrt(diag(m.SF2$vcov))
cbind(coef = m.SF2$coef, z = m.SF2$coef/sqrt(diag(m.SF2$vcov)))

}

}
 
 
 \source{
 \url{http://qed.econ.queensu.ca/jae/2014-v29.2/restrepo-tobon-kumbhakar/}.
 }
 
 \references{
Koetter, M., Kolari, J., and Spierdijk, L. (2012), Enjoying the quiet life under deregulation? Evidence from adjusted Lerner indices for U.S. banks. \emph{Review of Economics and Statistics}, \bold{94}, 2, 462--480.

Restrepo-Tobon, D. and Kumbhakar, S. (2014), Enjoying the quiet life under deregulation? Not Quite. \emph{Journal of Applied Econometrics}, \bold{29}, 2, 333--343.
 
 }
 

 
\keyword{datasets}
