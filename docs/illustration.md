# Illustration

## Data

> is a data frame containing selected variables for 500 U.S. commercial
> banks, randomly sampled from approximately 5000 banks, based on the
> dataset of Koetter et al. (2012) for year 2007. The dataset is
> provided solely for illustration and pedagogical purposes and is not
> suitable for empirical research.

``` r

  library(snreg)
  library(tidyverse)
#R>  ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
#R>  ✔ dplyr     1.1.4     ✔ readr     2.1.6
#R>  ✔ forcats   1.0.1     ✔ stringr   1.6.0
#R>  ✔ ggplot2   4.0.1     ✔ tibble    3.3.0
#R>  ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
#R>  ✔ purrr     1.2.0     
#R>  ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#R>  ✖ dplyr::filter() masks stats::filter()
#R>  ✖ dplyr::lag()    masks stats::lag()
#R>  ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
  data(banks07, package = "snreg")
  head(banks07)
#R>        id year        TC       Y1        Y2       W1       W2       W3        ER
#R>  1 152440 2007  4659.505  8661.84  58393.04 38.55422 37.75000 4.146092 0.1117789
#R>  2 544335 2007 10848.960 19492.91 146789.00 30.66785 45.13934 2.263886 0.1139121
#R>  3 548203 2007 14523.650 18630.49 219705.70 27.38797 56.02857 4.312672 0.1342984
#R>  4 568238 2007  1644.808 11902.50  18675.68 32.84768 50.47368 1.541316 0.1029891
#R>  5 651158 2007  3054.240 29322.21  38916.14 42.20033 52.23529 3.192056 0.1421004
#R>  6 705444 2007  6168.736 36568.03  67179.16 11.95771 41.46667 3.181241 0.1449880
#R>           LA        TA    ZSCORE  ZSCORE3     SDROA       LLP lnzscore lnzscore3
#R>  1 0.7781611  75039.78  29.49198 11.93898 0.4532092  732.4904 3.384118  2.479809
#R>  2 0.7679300 191148.92  62.97659 20.66493 0.2298561  305.9889 4.142763  3.028438
#R>  3 0.8880010 247416.05  38.99969 30.88642 0.3874364 1534.6521 3.663554  3.430317
#R>  4 0.5050412  36978.53  30.68124 24.53655 0.4063778    0.0000 3.423651  3.200164
#R>  5 0.5264002  73928.81  30.02672 28.42257 0.5496323    0.0000 3.402088  3.347184
#R>  6 0.5222276 128639.62 110.67650 37.66002 0.1502981    0.0000 4.706612  3.628599
#R>       lnsdroa     scope ms_county
#R>  1 -0.7914014 0.4801792 0.3466652
#R>  2 -1.4703018 0.5254947 0.8830606
#R>  3 -0.9482036 0.8897177 1.1430008
#R>  4 -0.9004720 0.3917692 0.1708316
#R>  5 -0.5985058 0.5523984 0.3415328
#R>  6 -1.8951346 0.4073012 0.5942832
```

## Specification

> Define the specification (formula) that will be used:

``` r

# Translog cost function specification
spe.tl <- log(TC) ~ (log(Y1) + log(Y2) + log(W1) + log(W2))^2 +
  I(0.5 * log(Y1)^2) + I(0.5 * log(Y2)^2) +
  I(0.5 * log(W1)^2) + I(0.5 * log(W2)^2)
```

## Linear Regression via MLE

To estimate simple OLS using MLE

``` r

# -------------------------------------------------------------
# Specification 1: homoskedastic noise (ln.var.v = NULL)
# -------------------------------------------------------------
formSV <- NULL

m1 <- lm.mle(
  formula   = spe.tl,
  data      = banks07,
  ln.var.v  = formSV
)
#R>   theta0:
#R>           (Intercept)              log(Y1)              log(Y2) 
#R>          -4.386991731          0.415874408          0.370050861 
#R>               log(W1)              log(W2)   I(0.5 * log(Y1)^2) 
#R>           0.524152369          1.342097991          0.047454756 
#R>    I(0.5 * log(Y2)^2)   I(0.5 * log(W1)^2)   I(0.5 * log(W2)^2) 
#R>           0.077484315          0.034996665         -0.226281101 
#R>       log(Y1):log(Y2)      log(Y1):log(W1)      log(Y1):log(W2) 
#R>          -0.057609631         -0.020462972         -0.005676474 
#R>       log(Y2):log(W1)      log(Y2):log(W2)      log(W1):log(W2) 
#R>           0.013639077          0.018840745         -0.154925076 
#R>  lnVARv0i_(Intercept) 
#R>          -3.687587268 
#R>  
#R>  -------------- Regression with normal errors: ---------------
#R>  
#R>  initial  value -212.427550 
#R>  iter   1 value -212.427550
#R>  final  value -212.427550 
#R>  converged
#R>                         Estimate    Std.Err  Z value    Pr(>z)    
#R>  (Intercept)          -4.3869917  3.0631686  -1.4322 0.1520939    
#R>  log(Y1)               0.4158744  0.1494204   2.7832 0.0053817 ** 
#R>  log(Y2)               0.3700509  0.2827874   1.3086 0.1906756    
#R>  log(W1)               0.5241524  0.2828127   1.8534 0.0638314 .  
#R>  log(W2)               1.3420980  1.2504546   1.0733 0.2831419    
#R>  I(0.5 * log(Y1)^2)    0.0474548  0.0060002   7.9088 2.665e-15 ***
#R>  I(0.5 * log(Y2)^2)    0.0774843  0.0225416   3.4374 0.0005873 ***
#R>  I(0.5 * log(W1)^2)    0.0349967  0.0273035   1.2818 0.1999249    
#R>  I(0.5 * log(W2)^2)   -0.2262811  0.3426730  -0.6603 0.5090349    
#R>  log(Y1):log(Y2)      -0.0576096  0.0113060  -5.0955 3.479e-07 ***
#R>  log(Y1):log(W1)      -0.0204630  0.0121495  -1.6843 0.0921316 .  
#R>  log(Y1):log(W2)      -0.0056765  0.0386938  -0.1467 0.8833669    
#R>  log(Y2):log(W1)       0.0136391  0.0161416   0.8450 0.3981301    
#R>  log(Y2):log(W2)       0.0188407  0.0596386   0.3159 0.7520666    
#R>  log(W1):log(W2)      -0.1549251  0.0695098  -2.2288 0.0258257 *  
#R>  lnVARv0i_(Intercept) -3.6875873  0.0632455 -58.3059 < 2.2e-16 ***
#R>  ---
#R>  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#R>  _____________________________________________________________

coef(m1)
#R>           (Intercept)              log(Y1)              log(Y2) 
#R>          -4.386991731          0.415874408          0.370050861 
#R>               log(W1)              log(W2)   I(0.5 * log(Y1)^2) 
#R>           0.524152369          1.342097991          0.047454756 
#R>    I(0.5 * log(Y2)^2)   I(0.5 * log(W1)^2)   I(0.5 * log(W2)^2) 
#R>           0.077484315          0.034996665         -0.226281101 
#R>       log(Y1):log(Y2)      log(Y1):log(W1)      log(Y1):log(W2) 
#R>          -0.057609631         -0.020462972         -0.005676474 
#R>       log(Y2):log(W1)      log(Y2):log(W2)      log(W1):log(W2) 
#R>           0.013639077          0.018840745         -0.154925076 
#R>  lnVARv0i_(Intercept) 
#R>          -3.687587268


# -------------------------------------------------------------
# Specification 2: heteroskedastic noise (variance depends on TA)
# -------------------------------------------------------------
formSV <- ~ log(TA)

m2 <- lm.mle(
  formula   = spe.tl,
  data      = banks07,
  ln.var.v  = formSV
)
#R>   theta0:
#R>           (Intercept)              log(Y1)              log(Y2) 
#R>          -4.386991731          0.415874408          0.370050861 
#R>               log(W1)              log(W2)   I(0.5 * log(Y1)^2) 
#R>           0.524152369          1.342097991          0.047454756 
#R>    I(0.5 * log(Y2)^2)   I(0.5 * log(W1)^2)   I(0.5 * log(W2)^2) 
#R>           0.077484315          0.034996665         -0.226281101 
#R>       log(Y1):log(Y2)      log(Y1):log(W1)      log(Y1):log(W2) 
#R>          -0.057609631         -0.020462972         -0.005676474 
#R>       log(Y2):log(W1)      log(Y2):log(W2)      log(W1):log(W2) 
#R>           0.013639077          0.018840745         -0.154925076 
#R>  lnVARv0i_(Intercept)     lnVARv0i_log(TA) 
#R>          -3.687587268          0.000000000 
#R>  
#R>  -------------- Regression with normal errors: ---------------
#R>  
#R>  initial  value -212.427550 
#R>  iter   2 value -212.435637
#R>  iter   2 value -212.435637
#R>  iter   2 value -212.435638
#R>  final  value -212.435638 
#R>  converged
#R>                          Estimate     Std.Err Z value    Pr(>z)    
#R>  (Intercept)          -4.38699173  3.07286157 -1.4277 0.1533907    
#R>  log(Y1)               0.41587442  0.14997873  2.7729 0.0055561 ** 
#R>  log(Y2)               0.37005087  0.29216784  1.2666 0.2053093    
#R>  log(W1)               0.52415237  0.28241172  1.8560 0.0634555 .  
#R>  log(W2)               1.34209799  1.25258014  1.0715 0.2839596    
#R>  I(0.5 * log(Y1)^2)    0.04745481  0.00604526  7.8499 4.219e-15 ***
#R>  I(0.5 * log(Y2)^2)    0.07748436  0.02254101  3.4375 0.0005871 ***
#R>  I(0.5 * log(W1)^2)    0.03499667  0.02724489  1.2845 0.1989592    
#R>  I(0.5 * log(W2)^2)   -0.22628109  0.34575683 -0.6545 0.5128209    
#R>  log(Y1):log(Y2)      -0.05760953  0.01151446 -5.0032 5.638e-07 ***
#R>  log(Y1):log(W1)      -0.02046294  0.01213895 -1.6857 0.0918485 .  
#R>  log(Y1):log(W2)      -0.00567644  0.03861266 -0.1470 0.8831243    
#R>  log(Y2):log(W1)       0.01363911  0.01617084  0.8434 0.3989832    
#R>  log(Y2):log(W2)       0.01884078  0.06006921  0.3137 0.7537860    
#R>  log(W1):log(W2)      -0.15492507  0.06936169 -2.2336 0.0255105 *  
#R>  lnVARv0i_(Intercept) -3.68758724  1.00965557 -3.6523 0.0002599 ***
#R>  lnVARv0i_log(TA)     -0.00036448  0.08630781 -0.0042 0.9966305    
#R>  ---
#R>  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#R>  _____________________________________________________________

coef(m2)
#R>           (Intercept)              log(Y1)              log(Y2) 
#R>         -4.3869917305         0.4158744170         0.3700508703 
#R>               log(W1)              log(W2)   I(0.5 * log(Y1)^2) 
#R>          0.5241523719         1.3420979943         0.0474548059 
#R>    I(0.5 * log(Y2)^2)   I(0.5 * log(W1)^2)   I(0.5 * log(W2)^2) 
#R>          0.0774843635         0.0349966695        -0.2262810941 
#R>       log(Y1):log(Y2)      log(Y1):log(W1)      log(Y1):log(W2) 
#R>         -0.0576095323        -0.0204629424        -0.0056764377 
#R>       log(Y2):log(W1)      log(Y2):log(W2)      log(W1):log(W2) 
#R>          0.0136391061         0.0188407818        -0.1549250651 
#R>  lnVARv0i_(Intercept)     lnVARv0i_log(TA) 
#R>         -3.6875872440        -0.0003644791
```

## Linear Regression with Skew-Normal Errors

> `snreg` fits a linear regression model where the disturbance term
> follows a skew-normal distribution.

``` r

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
#R>   theta0:
#R>           (Intercept)              log(Y1)              log(Y2) 
#R>          -4.386991731          0.415874408          0.370050861 
#R>               log(W1)              log(W2)   I(0.5 * log(Y1)^2) 
#R>           0.524152369          1.342097991          0.047454756 
#R>    I(0.5 * log(Y2)^2)   I(0.5 * log(W1)^2)   I(0.5 * log(W2)^2) 
#R>           0.077484315          0.034996665         -0.226281101 
#R>       log(Y1):log(Y2)      log(Y1):log(W1)      log(Y1):log(W2) 
#R>          -0.057609631         -0.020462972         -0.005676474 
#R>       log(Y2):log(W1)      log(Y2):log(W2)      log(W1):log(W2) 
#R>           0.013639077          0.018840745         -0.154925076 
#R>  lnVARv0i_(Intercept) Skew_v0i_(Intercept) 
#R>          -3.687587268          3.000000000 
#R>  
#R>  -------------- Regression with skewed errors: ---------------
#R>  
#R>  Iteration   0 (at starting values):       log likelihood = 5.146960520337
#R>  
#R>  Iteration   1 (hessian is provided,   1 in total):   log likelihood = 201.5619797261
#R>  Iteration   2 (hessian is provided,   2 in total):   log likelihood = 218.5157385369
#R>  Iteration   3 (hessian is provided,   3 in total):   log likelihood = 219.8304053615
#R>  Iteration   4 (hessian is provided,   4 in total):   log likelihood = 221.1570395072
#R>  Iteration   5 (hessian is provided,   5 in total):   log likelihood = 221.4173977326
#R>  Iteration   6 (hessian is provided,   6 in total):   log likelihood = 221.45910305
#R>  Iteration   7 (hessian is provided,   7 in total):   log likelihood = 221.4742381631
#R>  Iteration   8 (hessian is provided,   8 in total):   log likelihood = 221.4742892007
#R>  Iteration   9 (hessian is provided,   9 in total):   log likelihood = 221.4772909678
#R>  
#R>  Convergence given g inv(H) g' = 1.147956e-05 < lmtol
#R>  
#R>  Final log likelihood = 221.4772909678
#R>  
#R>                          Estimate     Std.Err  Z value    Pr(>z)    
#R>  (Intercept)          -5.88439482  3.00884171  -1.9557  0.050500 .  
#R>  log(Y1)               0.43986913  0.15724371   2.7974  0.005152 ** 
#R>  log(Y2)               0.50497847  0.27574767   1.8313  0.067055 .  
#R>  log(W1)               0.54680112  0.27757332   1.9699  0.048846 *  
#R>  log(W2)               1.62221319  1.19336826   1.3594  0.174034    
#R>  I(0.5 * log(Y1)^2)    0.04152960  0.00579831   7.1624 7.929e-13 ***
#R>  I(0.5 * log(Y2)^2)    0.06775246  0.02251267   3.0095  0.002617 ** 
#R>  I(0.5 * log(W1)^2)    0.03985513  0.02624261   1.5187  0.128833    
#R>  I(0.5 * log(W2)^2)   -0.28421570  0.32037935  -0.8871  0.375013    
#R>  log(Y1):log(Y2)      -0.05595449  0.01162258  -4.8143 1.477e-06 ***
#R>  log(Y1):log(W1)      -0.02195339  0.01251279  -1.7545  0.079349 .  
#R>  log(Y1):log(W2)       0.00032291  0.04106930   0.0079  0.993727    
#R>  log(Y2):log(W1)       0.01294569  0.01577160   0.8208  0.411747    
#R>  log(Y2):log(W2)       0.00975499  0.05798118   0.1682  0.866391    
#R>  log(W1):log(W2)      -0.15803381  0.06620501  -2.3870  0.016985 *  
#R>  lnVARv0i_(Intercept) -3.01076847  0.11366937 -26.4871 < 2.2e-16 ***
#R>  Skew_v0i_(Intercept)  1.86603587  0.32264117   5.7836 7.311e-09 ***
#R>  ---
#R>  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#R>  _____________________________________________________________

coef(m1)
#R>           (Intercept)              log(Y1)              log(Y2) 
#R>          -5.884394817          0.439869125          0.504978469 
#R>               log(W1)              log(W2)   I(0.5 * log(Y1)^2) 
#R>           0.546801118          1.622213193          0.041529597 
#R>    I(0.5 * log(Y2)^2)   I(0.5 * log(W1)^2)   I(0.5 * log(W2)^2) 
#R>           0.067752458          0.039855132         -0.284215702 
#R>       log(Y1):log(Y2)      log(Y1):log(W1)      log(Y1):log(W2) 
#R>          -0.055954492         -0.021953387          0.000322907 
#R>       log(Y2):log(W1)      log(Y2):log(W2)      log(W1):log(W2) 
#R>           0.012945689          0.009754988         -0.158033814 
#R>  lnVARv0i_(Intercept) Skew_v0i_(Intercept) 
#R>          -3.010768472          1.866035868


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
#R>   theta0:
#R>                (Intercept)                   log(Y1)                   log(Y2) 
#R>               -4.386991731               0.415874408               0.370050861 
#R>                    log(W1)                   log(W2)        I(0.5 * log(Y1)^2) 
#R>                0.524152369               1.342097991               0.047454756 
#R>         I(0.5 * log(Y2)^2)        I(0.5 * log(W1)^2)        I(0.5 * log(W2)^2) 
#R>                0.077484315               0.034996665              -0.226281101 
#R>            log(Y1):log(Y2)           log(Y1):log(W1)           log(Y1):log(W2) 
#R>               -0.057609631              -0.020462972              -0.005676474 
#R>            log(Y2):log(W1)           log(Y2):log(W2)           log(W1):log(W2) 
#R>                0.013639077               0.018840745              -0.154925076 
#R>       lnVARv0i_(Intercept) lnVARv0i_lnVARv0i_log(TA)      Skew_v0i_(Intercept) 
#R>               -3.687587268               0.000000000               3.000000000 
#R>       Skew_v0i_Skew_v0i_ER 
#R>                0.000000000 
#R>  
#R>  ----------------- Regression with skewed errors: -----------------
#R>  
#R>  Iteration   0 (at starting values):       log likelihood = 5.146960520337
#R>  
#R>  Iteration   1 (hessian is provided,   1 in total):   log likelihood = 204.6988401214
#R>  Iteration   2 (hessian is provided,   2 in total):   log likelihood = 221.6390613234
#R>  Iteration   3 (hessian is provided,   3 in total):   log likelihood = 224.8771023643
#R>  Iteration   4 (hessian is provided,   4 in total):   log likelihood = 225.6867346145
#R>  Iteration   5 (hessian is provided,   5 in total):   log likelihood = 225.9380532273
#R>  Iteration   6 (hessian is provided,   6 in total):   log likelihood = 226.0484166119
#R>  Iteration   7 (hessian is provided,   7 in total):   log likelihood = 226.0721225638
#R>  Iteration   8 (hessian is provided,   8 in total):   log likelihood = 226.0820635056
#R>  Iteration   9 (hessian is provided,   9 in total):   log likelihood = 226.0835072457
#R>  Iteration  10 (hessian is provided,  10 in total):   log likelihood = 226.0909506421
#R>  
#R>  Convergence given g inv(H) g' = 9.357213e-05 < lmtol
#R>  
#R>  Final log likelihood = 226.0909506421
#R>  
#R>                               Estimate     Std.Err Z value    Pr(>z)    
#R>  (Intercept)               -5.71501334  2.98991742 -1.9114  0.055950 .  
#R>  log(Y1)                    0.45216186  0.15708195  2.8785  0.003996 ** 
#R>  log(Y2)                    0.64671047  0.28153250  2.2971  0.021613 *  
#R>  log(W1)                    0.48917765  0.27664356  1.7683  0.077018 .  
#R>  log(W2)                    1.15530138  1.17400483  0.9841  0.325082    
#R>  I(0.5 * log(Y1)^2)         0.04070906  0.00590546  6.8935 5.445e-12 ***
#R>  I(0.5 * log(Y2)^2)         0.06170967  0.02265802  2.7235  0.006459 ** 
#R>  I(0.5 * log(W1)^2)         0.04271441  0.02667730  1.6012  0.109343    
#R>  I(0.5 * log(W2)^2)        -0.15999322  0.31994376 -0.5001  0.617028    
#R>  log(Y1):log(Y2)           -0.05757860  0.01197968 -4.8064 1.537e-06 ***
#R>  log(Y1):log(W1)           -0.01887894  0.01276125 -1.4794  0.139035    
#R>  log(Y1):log(W2)            0.00165155  0.04117676  0.0401  0.968006    
#R>  log(Y2):log(W1)            0.00516054  0.01583258  0.3259  0.744467    
#R>  log(Y2):log(W2)           -0.00022261  0.05922730 -0.0038  0.997001    
#R>  log(W1):log(W2)           -0.13107647  0.06626509 -1.9781  0.047922 *  
#R>  lnVARv0i_(Intercept)      -0.89452888  1.03829604 -0.8615  0.388943    
#R>  lnVARv0i_lnVARv0i_log(TA) -0.17730037  0.08736480 -2.0294  0.042415 *  
#R>  Skew_v0i_(Intercept)       3.10567750  0.76059557  4.0832 4.442e-05 ***
#R>  Skew_v0i_Skew_v0i_ER      -9.55208958  4.87607389 -1.9590  0.050116 .  
#R>  ---
#R>  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#R>  __________________________________________________________________

coef(m2)
#R>                (Intercept)                   log(Y1)                   log(Y2) 
#R>              -5.7150133374              0.4521618623              0.6467104718 
#R>                    log(W1)                   log(W2)        I(0.5 * log(Y1)^2) 
#R>               0.4891776489              1.1553013814              0.0407090640 
#R>         I(0.5 * log(Y2)^2)        I(0.5 * log(W1)^2)        I(0.5 * log(W2)^2) 
#R>               0.0617096653              0.0427144069             -0.1599932243 
#R>            log(Y1):log(Y2)           log(Y1):log(W1)           log(Y1):log(W2) 
#R>              -0.0575785975             -0.0188789363              0.0016515541 
#R>            log(Y2):log(W1)           log(Y2):log(W2)           log(W1):log(W2) 
#R>               0.0051605363             -0.0002226092             -0.1310764688 
#R>       lnVARv0i_(Intercept) lnVARv0i_lnVARv0i_log(TA)      Skew_v0i_(Intercept) 
#R>              -0.8945288842             -0.1773003673              3.1056774984 
#R>       Skew_v0i_Skew_v0i_ER 
#R>              -9.5520895773
```

## Stochastic Frontier Model with a Skew-Normally Distributed Error Term

> `snsf` performs maximum likelihood estimation of the parameters and
> technical or cost efficiencies in a Stochastic Frontier Model with a
> skew-normally distributed error term.

``` r

myprod <- FALSE

# -------------------------------------------------------------
# Specification 1: homoskedastic & symmetric
# -------------------------------------------------------------
formSV <- NULL   # variance equation
formSK <- NULL   # skewness equation
formSU <- NULL   # inefficiency equation (unused here)

m1 <- snsf(
  formula  = spe.tl,
  data     = banks07,
  prod     = myprod,
  ln.var.v = formSV,
  skew.v   = formSK
)
#R>  
#R>  -------------- Regression with skewed errors: ---------------
#R>  
#R>  initial  value -105.392823 
#R>  iter   2 value -105.561746
#R>  iter   3 value -105.568000
#R>  iter   4 value -105.836407
#R>  iter   5 value -106.255928
#R>  iter   6 value -107.989401
#R>  iter   7 value -129.390401
#R>  iter   8 value -148.614164
#R>  iter   8 value -148.614164
#R>  iter   9 value -148.627674
#R>  iter  10 value -148.628304
#R>  iter  10 value -148.628304
#R>  iter  10 value -148.628304
#R>  final  value -148.628304 
#R>  converged
#R>                          Estimate     Std.Err  Z value    Pr(>z)    
#R>  X(Intercept)         -4.40500089  3.80254827  -1.1584  0.246687    
#R>  Xlog(Y1)              0.45706967  0.18307828   2.4966  0.012540 *  
#R>  Xlog(Y2)              0.31560584  0.33532057   0.9412  0.346599    
#R>  Xlog(W1)              0.48827077  0.38637994   1.2637  0.206335    
#R>  Xlog(W2)              1.29215252  1.70829901   0.7564  0.449411    
#R>  XI(0.5 * log(Y1)^2)   0.06252731  0.00718640   8.7008 < 2.2e-16 ***
#R>  XI(0.5 * log(Y2)^2)   0.19680294  0.03130033   6.2876 3.225e-10 ***
#R>  XI(0.5 * log(W1)^2)   0.05041646  0.03723206   1.3541  0.175700    
#R>  XI(0.5 * log(W2)^2)  -0.28522159  0.46697561  -0.6108  0.541342    
#R>  Xlog(Y1):log(Y2)     -0.14258698  0.00610987 -23.3372 < 2.2e-16 ***
#R>  Xlog(Y1):log(W1)      0.00029592  0.01705705   0.0173  0.986158    
#R>  Xlog(Y1):log(W2)      0.17610934  0.05558196   3.1685  0.001532 ** 
#R>  Xlog(Y2):log(W1)      0.01454480  0.02162340   0.6726  0.501175    
#R>  Xlog(Y2):log(W2)     -0.09354141  0.07710129  -1.2132  0.225043    
#R>  Xlog(W1):log(W2)     -0.21197450  0.09475124  -2.2372  0.025275 *  
#R>  lnVARv0i_(Intercept) -2.66540968  0.11013951 -24.2003 < 2.2e-16 ***
#R>  Skew_v0i_(Intercept)  0.98854923  0.03434255  28.7850 < 2.2e-16 ***
#R>  sv                    0.26376286  0.00766250  34.4226 < 2.2e-16 ***
#R>  ---
#R>  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#R>  _____________________________________________________________
#R>   theta0:
#R>           (Intercept)              log(Y1)              log(Y2) 
#R>         -4.4050008900         0.4570696682         0.3156058426 
#R>               log(W1)              log(W2)   I(0.5 * log(Y1)^2) 
#R>          0.4882707672         1.2921525200         0.0625273146 
#R>    I(0.5 * log(Y2)^2)   I(0.5 * log(W1)^2)   I(0.5 * log(W2)^2) 
#R>          0.1968029394         0.0504164557        -0.2852215911 
#R>       log(Y1):log(Y2)      log(Y1):log(W1)      log(Y1):log(W2) 
#R>         -0.1425869800         0.0002959179         0.1761093362 
#R>       log(Y2):log(W1)      log(Y2):log(W2)      log(W1):log(W2) 
#R>          0.0145448021        -0.0935414140        -0.2119744969 
#R>  lnVARv0i_(Intercept) Skew_v0i_(Intercept) lnVARu0i_(Intercept) 
#R>         -4.2020745955         0.3982617494        -4.5984104618 
#R>  
#R>  ---------------------- The main model: ----------------------
#R>  
#R>  Iteration   0 (at starting values):       log likelihood = 21.64346032248
#R>  
#R>  Iteration   1 (hessian is provided,   1 in total):   log likelihood = 186.4880063145
#R>  Iteration   2 (hessian is provided,   2 in total):   log likelihood = 214.3430643968
#R>  Iteration   3 (hessian is provided,   3 in total):   log likelihood = 220.957921695
#R>  Iteration   4 (hessian is provided,   4 in total):   log likelihood = 223.1040419712
#R>  Iteration   5 (hessian is provided,   5 in total):   log likelihood = 224.3725358921
#R>  Iteration   6 (hessian is provided,   6 in total):   log likelihood = 224.6133195507
#R>  Iteration   7 (hessian is provided,   7 in total):   log likelihood = 224.7096271178
#R>  Iteration   8 (hessian is provided,   8 in total):   log likelihood = 224.7298349674
#R>  Iteration   9 (hessian is provided,   9 in total):   log likelihood = 224.8005650898
#R>  Iteration  10 (hessian is provided,  10 in total):   log likelihood = 224.8011109934
#R>  Iteration  11 (hessian is provided,  11 in total):   log likelihood = 224.8019248405
#R>  Iteration  12 (hessian is provided,  12 in total):   log likelihood = 224.8021079368
#R>  Iteration  13 (hessian is provided,  13 in total):   log likelihood = 224.8054321542
#R>  Iteration  14 (hessian is provided,  14 in total):   log likelihood = 224.8420044777
#R>  Iteration  15 (hessian is provided,  15 in total):   log likelihood = 225.4135588864
#R>  Iteration  16 (hessian is provided,  16 in total):   log likelihood = 225.5282232035
#R>  Iteration  17 (hessian is provided,  17 in total):   log likelihood = 225.533351752
#R>  Iteration  18 (hessian is provided,  18 in total):   log likelihood = 225.5368333384
#R>  Iteration  19 (hessian is provided,  19 in total):   log likelihood = 225.5371105355
#R>  Iteration  20 (hessian is provided,  20 in total):   log likelihood = 225.539300346
#R>  Iteration  21 (hessian is provided,  21 in total):   log likelihood = 225.539317247
#R>  Iteration  22 (hessian is provided,  22 in total):   log likelihood = 225.539324279
#R>  
#R>  Convergence given g inv(H) g' = 1.027857e-06 < lmtol
#R>  
#R>  Final log likelihood = 225.539324279
#R>  
#R>  
#R>  Log likelihood maximization completed in
#R>  0.18 seconds
#R>  _____________________________________________________________
#R>  
#R>  Cross-sectional stochastic (cost) frontier model
#R>  
#R>  Distributional assumptions
#R>  
#R>    Component      Distribution Assumption   
#R>  1 Random noise:  skew normal  homoskedastic
#R>  2 Inefficiency:  exponential  homoskedastic
#R>  
#R>  Number of observations = 500
#R>  
#R>  -------------------- Estimation results: --------------------
#R>  
#R>                            Coef.        SE       z       P>|z|
#R>  _____________________________________________________________
#R>  Frontier
#R>                                                                   
#R>  (Intercept)             -6.3606     3.0602   -2.08     0.0377 *  
#R>  log(Y1)                  0.4425     0.1767    2.50     0.0123 *  
#R>  log(Y2)                  0.5660     0.2753    2.06     0.0398 *  
#R>  log(W1)                  0.4983     0.2835    1.76     0.0788 .  
#R>  log(W2)                  1.6557     1.1889    1.39     0.1637    
#R>  I(0.5 * log(Y1)^2)       0.0449     0.0066    6.76     0.0000 ***
#R>  I(0.5 * log(Y2)^2)       0.0650     0.0222    2.93     0.0034 ** 
#R>  I(0.5 * log(W1)^2)       0.0363     0.0258    1.41     0.1597    
#R>  I(0.5 * log(W2)^2)      -0.3099     0.3170   -0.98     0.3283    
#R>  log(Y1):log(Y2)         -0.0586     0.0117   -4.99     0.0000 ***
#R>  log(Y1):log(W1)         -0.0227     0.0127   -1.78     0.0750 .  
#R>  log(Y1):log(W2)         -0.0006     0.0420   -0.01     0.9895    
#R>  log(Y2):log(W1)          0.0117     0.0152    0.77     0.4393    
#R>  log(Y2):log(W2)          0.0106     0.0563    0.19     0.8505    
#R>  log(W1):log(W2)         -0.1375     0.0682   -2.02     0.0438 *  
#R>  -------------------------------------------------------------
#R>  Variance of the random noise component: log(sigma_v^2)
#R>                                                                   
#R>  lnVARv0i_(Intercept)    -3.8207     0.2255  -16.94     0.0000 ***
#R>  -------------------------------------------------------------
#R>  Skewness of the random noise component: `alpha`
#R>                                                                   
#R>  Skew_v0i_(Intercept)    -1.5007     0.8755   -1.71     0.0865 .  
#R>  -------------------------------------------------------------
#R>  Inefficiency component: log(sigma_u^2)
#R>                                                                   
#R>  lnVARu0i_(Intercept)    -4.3422     0.2398  -18.11     0.0000 ***
#R>  _____________________________________________________________
#R>  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#R>  
#R>  --------------- Summary of cost efficiencies: ---------------
#R>  
#R>                  JLMS:= exp(-E[ui|ei])
#R>  
#R>        TE_JLMS
#R>  Obs       500
#R>  NAs         0
#R>  Mean   0.8955
#R>  StDev  0.0715
#R>  IQR    0.0645
#R>  Min    0.4559
#R>  5%     0.7420
#R>  10%    0.7959
#R>  25%    0.8780
#R>  50%    0.9208
#R>  75%    0.9425
#R>  90%    0.9520
#R>  Max    0.9691

coef(m1)
#R>           (Intercept)              log(Y1)              log(Y2) 
#R>         -6.3606134654         0.4425312174         0.5660453675 
#R>               log(W1)              log(W2)   I(0.5 * log(Y1)^2) 
#R>          0.4983319531         1.6557355502         0.0448519992 
#R>    I(0.5 * log(Y2)^2)   I(0.5 * log(W1)^2)   I(0.5 * log(W2)^2) 
#R>          0.0649995029         0.0362747554        -0.3099280819 
#R>       log(Y1):log(Y2)      log(Y1):log(W1)      log(Y1):log(W2) 
#R>         -0.0586056701        -0.0226917826        -0.0005509121 
#R>       log(Y2):log(W1)      log(Y2):log(W2)      log(W1):log(W2) 
#R>          0.0117430476         0.0106095399        -0.1375404433 
#R>  lnVARv0i_(Intercept) Skew_v0i_(Intercept) lnVARu0i_(Intercept) 
#R>         -3.8206601915        -1.5006628040        -4.3422317628


# -------------------------------------------------------------
# Specification 2: heteroskedastic + skewed noise
# -------------------------------------------------------------
formSV <- ~ log(TA)      # heteroskedastic variance
formSK <- ~ ER           # skewness driver
formSU <- ~ LA + ER      # inefficiency

m2 <- snsf(
  formula  = spe.tl,
  data     = banks07,
  prod     = myprod,
  ln.var.v = formSV,
  skew.v   = formSK
)
#R>  
#R>  -------------- Regression with skewed errors: ---------------
#R>  
#R>  initial  value -105.392823 
#R>  iter   2 value -105.561746
#R>  iter   3 value -105.568000
#R>  iter   4 value -105.836407
#R>  iter   5 value -106.255928
#R>  iter   6 value -107.989401
#R>  iter   7 value -129.390401
#R>  iter   8 value -148.614164
#R>  iter   8 value -148.614164
#R>  iter   9 value -148.627674
#R>  iter  10 value -148.628304
#R>  iter  10 value -148.628304
#R>  iter  10 value -148.628304
#R>  final  value -148.628304 
#R>  converged
#R>                          Estimate     Std.Err  Z value    Pr(>z)    
#R>  X(Intercept)         -4.40500089  3.80254827  -1.1584  0.246687    
#R>  Xlog(Y1)              0.45706967  0.18307828   2.4966  0.012540 *  
#R>  Xlog(Y2)              0.31560584  0.33532057   0.9412  0.346599    
#R>  Xlog(W1)              0.48827077  0.38637994   1.2637  0.206335    
#R>  Xlog(W2)              1.29215252  1.70829901   0.7564  0.449411    
#R>  XI(0.5 * log(Y1)^2)   0.06252731  0.00718640   8.7008 < 2.2e-16 ***
#R>  XI(0.5 * log(Y2)^2)   0.19680294  0.03130033   6.2876 3.225e-10 ***
#R>  XI(0.5 * log(W1)^2)   0.05041646  0.03723206   1.3541  0.175700    
#R>  XI(0.5 * log(W2)^2)  -0.28522159  0.46697561  -0.6108  0.541342    
#R>  Xlog(Y1):log(Y2)     -0.14258698  0.00610987 -23.3372 < 2.2e-16 ***
#R>  Xlog(Y1):log(W1)      0.00029592  0.01705705   0.0173  0.986158    
#R>  Xlog(Y1):log(W2)      0.17610934  0.05558196   3.1685  0.001532 ** 
#R>  Xlog(Y2):log(W1)      0.01454480  0.02162340   0.6726  0.501175    
#R>  Xlog(Y2):log(W2)     -0.09354141  0.07710129  -1.2132  0.225043    
#R>  Xlog(W1):log(W2)     -0.21197450  0.09475124  -2.2372  0.025275 *  
#R>  lnVARv0i_(Intercept) -2.66540968  0.11013951 -24.2003 < 2.2e-16 ***
#R>  Skew_v0i_(Intercept)  0.98854923  0.03434255  28.7850 < 2.2e-16 ***
#R>  sv                    0.26376286  0.00766250  34.4226 < 2.2e-16 ***
#R>  ---
#R>  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#R>  _____________________________________________________________
#R>   theta0:
#R>           (Intercept)              log(Y1)              log(Y2) 
#R>         -4.405001e+00         4.570697e-01         3.156058e-01 
#R>               log(W1)              log(W2)   I(0.5 * log(Y1)^2) 
#R>          4.882708e-01         1.292153e+00         6.252731e-02 
#R>    I(0.5 * log(Y2)^2)   I(0.5 * log(W1)^2)   I(0.5 * log(W2)^2) 
#R>          1.968029e-01         5.041646e-02        -2.852216e-01 
#R>       log(Y1):log(Y2)      log(Y1):log(W1)      log(Y1):log(W2) 
#R>         -1.425870e-01         2.959179e-04         1.761093e-01 
#R>       log(Y2):log(W1)      log(Y2):log(W2)      log(W1):log(W2) 
#R>          1.454480e-02        -9.354141e-02        -2.119745e-01 
#R>  lnVARv0i_(Intercept)     lnVARv0i_log(TA) Skew_v0i_(Intercept) 
#R>         -4.202075e+00         2.212626e-16         3.982617e-01 
#R>           Skew_v0i_ER lnVARu0i_(Intercept) 
#R>          0.000000e+00        -4.598410e+00 
#R>  
#R>  ---------------------- The main model: ----------------------
#R>  
#R>  Iteration   0 (at starting values):       log likelihood = 21.64346032248
#R>  
#R>  Iteration   1 (hessian is provided,   1 in total):   log likelihood = 133.2018982456
#R>  Iteration   2 (hessian is provided,   2 in total):   log likelihood = 183.2824024606
#R>  Iteration   3 (hessian is provided,   3 in total):   log likelihood = 214.7412845581
#R>  Iteration   4 (hessian is provided,   4 in total):   log likelihood = 223.1763227676
#R>  Iteration   5 (hessian is provided,   5 in total):   log likelihood = 226.655477596
#R>  Iteration   6 (hessian is provided,   6 in total):   log likelihood = 226.7105002202
#R>  Iteration   7 (hessian is provided,   7 in total):   log likelihood = 226.8074129256
#R>  Iteration   8 (hessian is provided,   8 in total):   log likelihood = 226.8082101988
#R>  Iteration   9 (hessian is provided,   9 in total):   log likelihood = 226.8082151929
#R>  Iteration  10 (hessian is provided,  10 in total):   log likelihood = 226.808224564
#R>  
#R>  Convergence given g inv(H) g' = 2.529443e-06 < lmtol
#R>  
#R>  Final log likelihood = 226.808224564
#R>  
#R>  
#R>  Log likelihood maximization completed in
#R>  0.103 seconds
#R>  _____________________________________________________________
#R>  
#R>  Cross-sectional stochastic (cost) frontier model
#R>  
#R>  Distributional assumptions
#R>  
#R>    Component      Distribution Assumption     
#R>  1 Random noise:  skew normal  heteroskedastic
#R>  2 Inefficiency:  exponential  homoskedastic  
#R>  
#R>  Number of observations = 500
#R>  
#R>  -------------------- Estimation results: --------------------
#R>  
#R>                            Coef.        SE       z       P>|z|
#R>  _____________________________________________________________
#R>  Frontier
#R>                                                                   
#R>  (Intercept)             -5.9441     3.0205   -1.97     0.0491 *  
#R>  log(Y1)                  0.4589     0.1647    2.79     0.0053 ** 
#R>  log(Y2)                  0.6135     0.2820    2.18     0.0296 *  
#R>  log(W1)                  0.5168     0.2747    1.88     0.0599 .  
#R>  log(W2)                  1.2707     1.1934    1.06     0.2870    
#R>  I(0.5 * log(Y1)^2)       0.0414     0.0061    6.81     0.0000 ***
#R>  I(0.5 * log(Y2)^2)       0.0614     0.0226    2.72     0.0065 ** 
#R>  I(0.5 * log(W1)^2)       0.0439     0.0262    1.68     0.0938 .  
#R>  I(0.5 * log(W2)^2)      -0.2078     0.3247   -0.64     0.5223    
#R>  log(Y1):log(Y2)         -0.0575     0.0120   -4.79     0.0000 ***
#R>  log(Y1):log(W1)         -0.0215     0.0128   -1.68     0.0927 .  
#R>  log(Y1):log(W2)          0.0000     0.0417    0.00     0.9997    
#R>  log(Y2):log(W1)          0.0068     0.0158    0.43     0.6682    
#R>  log(Y2):log(W2)          0.0086     0.0589    0.15     0.8836    
#R>  log(W1):log(W2)         -0.1368     0.0661   -2.07     0.0384 *  
#R>  -------------------------------------------------------------
#R>  Variance of the random noise component: log(sigma_v^2)
#R>                                                                   
#R>  lnVARv0i_(Intercept)    -1.8839     1.6595   -1.14     0.2563    
#R>  lnVARv0i_log(TA)        -0.1558     0.1328   -1.17     0.2407    
#R>  -------------------------------------------------------------
#R>  Skewness of the random noise component: `alpha`
#R>                                                                   
#R>  Skew_v0i_(Intercept)     2.0043     1.0468    1.91     0.0555 .  
#R>  Skew_v0i_ER             -7.2632     5.0581   -1.44     0.1510    
#R>  -------------------------------------------------------------
#R>  Inefficiency component: log(sigma_u^2)
#R>                                                                   
#R>  lnVARu0i_(Intercept)    -4.6407     0.3416  -13.58     0.0000 ***
#R>  _____________________________________________________________
#R>  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#R>  
#R>  --------------- Summary of cost efficiencies: ---------------
#R>  
#R>                  JLMS:= exp(-E[ui|ei])
#R>  
#R>        TE_JLMS
#R>  Obs       500
#R>  NAs         0
#R>  Mean   0.9084
#R>  StDev  0.0564
#R>  IQR    0.0547
#R>  Min    0.5206
#R>  5%     0.8096
#R>  10%    0.8379
#R>  25%    0.8903
#R>  50%    0.9241
#R>  75%    0.9450
#R>  90%    0.9568
#R>  Max    0.9758

coef(m2)
#R>           (Intercept)              log(Y1)              log(Y2) 
#R>         -5.944090e+00         4.589482e-01         6.134700e-01 
#R>               log(W1)              log(W2)   I(0.5 * log(Y1)^2) 
#R>          5.167646e-01         1.270694e+00         4.138294e-02 
#R>    I(0.5 * log(Y2)^2)   I(0.5 * log(W1)^2)   I(0.5 * log(W2)^2) 
#R>          6.141001e-02         4.385504e-02        -2.077538e-01 
#R>       log(Y1):log(Y2)      log(Y1):log(W1)      log(Y1):log(W2) 
#R>         -5.745504e-02        -2.151307e-02         1.351571e-05 
#R>       log(Y2):log(W1)      log(Y2):log(W2)      log(W1):log(W2) 
#R>          6.752251e-03         8.627646e-03        -1.368348e-01 
#R>  lnVARv0i_(Intercept)     lnVARv0i_log(TA) Skew_v0i_(Intercept) 
#R>         -1.883867e+00        -1.558322e-01         2.004266e+00 
#R>           Skew_v0i_ER lnVARu0i_(Intercept) 
#R>         -7.263169e+00        -4.640750e+00
```

## Additional Resources

To be added
