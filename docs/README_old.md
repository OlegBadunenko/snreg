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
