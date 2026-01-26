## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----setup--------------------------------------------------------------------
library(snreg)
library(tidyverse)

## ----data---------------------------------------------------------------------
data(banks07, package = "snreg")
head(banks07)

## ----formula------------------------------------------------------------------
spe.tl <- log(TC) ~ (log(Y1) + log(Y2) + log(W1) + log(W2))^2 +
  I(0.5 * log(Y1)^2) + I(0.5 * log(Y2)^2) +
  I(0.5 * log(W1)^2) + I(0.5 * log(W2)^2)

## ----lmmle1-------------------------------------------------------------------
formSV <- NULL

m1 <- lm.mle(
  formula   = spe.tl,
  data      = banks07,
  ln.var.v  = formSV
)

coef(m1)

## ----lmmle2-------------------------------------------------------------------
formSV <- ~ log(TA)

m2 <- lm.mle(
  formula   = spe.tl,
  data      = banks07,
  ln.var.v  = formSV
)

coef(m2)

## ----snreg1-------------------------------------------------------------------
formSV <- NULL     # variance equation
formSK <- NULL     # skewness equation

m1 <- snreg(
  formula  = spe.tl,
  data     = banks07,
  ln.var.v = formSV,
  skew.v   = formSK
)

coef(m1)

## ----snreg2-------------------------------------------------------------------
formSV <- ~ log(TA)   # heteroskedasticity in v
formSK <- ~ ER        # skewness driven by equity ratio

m2 <- snreg(
  formula  = spe.tl,
  data     = banks07,
  ln.var.v = formSV,
  skew.v   = formSK
)

coef(m2)

## ----snsf1--------------------------------------------------------------------
myprod <- FALSE

formSV <- NULL   # variance equation
formSK <- NULL   # skewness equation

m1 <- snsf(
  formula  = spe.tl,
  data     = banks07,
  prod     = myprod,
  ln.var.v = formSV,
  skew.v   = formSK
)

coef(m1)

## ----snsf2--------------------------------------------------------------------
formSV <- ~ log(TA)      # heteroskedastic variance
formSK <- ~ ER           # skewness driver

m2 <- snsf(
  formula  = spe.tl,
  data     = banks07,
  prod     = myprod,
  ln.var.v = formSV,
  skew.v   = formSK
)

coef(m2)

