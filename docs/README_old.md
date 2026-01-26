# 📦 Installing snreg R Package from GitHub

## The framework

Repository for the R package `snreg`.

R commands `snreg` and `snsf` estimate models with skew-normal errors
written and maintained by Oleg Badunenko
(<oleg.badunenko@brunel.ac.uk>). The details are discussed in

Badunenko, O., & Henderson, D. J. (2023). Production analysis with
asymmetric noise. Journal of Productivity Analysis, 61(1), 1–18. [DOI
![](reference/figures/doi.png)](https://doi.org/10.1007/s11123-023-00680-5)

### 1. Install the `devtools` package

> Install `devtools` from CRAN (if you haven’t already):

``` r

install.packages("devtools")
```

### 2. Install the package from GitHub

> Install the `snreg` package from GitHub. In this code, you are
> installing the `snreg` package created by `OlegBadunenko`.

``` r

library(devtools)
install_github("OlegBadunenko/snreg")
```

### 3. Load the installed package

``` r

library(snreg)
```

## 💡 Notes & Tips

- Works identically across **R**, **RStudio**, **Windows**, **Mac**, and
  **Linux**.
- Some GitHub packages may already be available in environments like
  `npsf`.
- If installation fails, common causes include missing build tools,
  incorrect repo names, or network restrictions.

# Illustration

## Data

> is a data frame containing selected variables for 500 U.S. commercial
> banks, randomly sampled from approximately 5000 banks, based on the
> dataset of Koetter et al. (2012) for year 2007. The dataset is
> provided solely for illustration and pedagogical purposes and is not
> suitable for empirical research.

`{r data, eval = TRUE} library(snreg) library(tidyverse) data(banks07, package = "snreg") head(banks07)`

## Specification

> Define the specification (formula) that will be used:

`{r formula, eval = TRUE} # Translog cost function specification spe.tl <- log(TC) ~ (log(Y1) + log(Y2) + log(W1) + log(W2))^2 + I(0.5 * log(Y1)^2) + I(0.5 * log(Y2)^2) + I(0.5 * log(W1)^2) + I(0.5 * log(W2)^2)`

## Linear Regression via MLE

To estimate simple OLS using MLE

\`\`\`{r lmmle, eval = TRUE} \# ————————————————————- \# Specification
1: homoskedastic noise (ln.var.v = NULL) \# ————————————————————- formSV
\<- NULL

m1 \<- lm.mle( formula = spe.tl, data = banks07, ln.var.v = formSV )

coef(m1)

# ————————————————————-

# Specification 2: heteroskedastic noise (variance depends on TA)

# ————————————————————-

formSV \<- ~ log(TA)

m2 \<- lm.mle( formula = spe.tl, data = banks07, ln.var.v = formSV )

coef(m2)


    ## Linear Regression with Skew-Normal Errors

    > `snreg` fits a linear regression model where the disturbance term follows a skew-normal distribution.

    ```{r snreg, eval=TRUE}
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

    coef(m1)


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

    coef(m2)

## Stochastic Frontier Model with a Skew-Normally Distributed Error Term

> `snsf` performs maximum likelihood estimation of the parameters and
> technical or cost efficiencies in a Stochastic Frontier Model with a
> skew-normally distributed error term.

\`\`\`{r snsf, eval=TRUE} myprod \<- FALSE

# ————————————————————-

# Specification 1: homoskedastic & symmetric

# ————————————————————-

formSV \<- NULL \# variance equation formSK \<- NULL \# skewness
equation formSU \<- NULL \# inefficiency equation (unused here)

m1 \<- snsf( formula = spe.tl, data = banks07, prod = myprod, ln.var.v =
formSV, skew.v = formSK )

coef(m1)

# ————————————————————-

# Specification 2: heteroskedastic + skewed noise

# ————————————————————-

formSV \<- ~ log(TA) \# heteroskedastic variance formSK \<- ~ ER \#
skewness driver formSU \<- ~ LA + ER \# inefficiency

m2 \<- snsf( formula = spe.tl, data = banks07, prod = myprod, ln.var.v =
formSV, skew.v = formSK )

coef(m2) \`\`\`

## Additional Resources

To be added
