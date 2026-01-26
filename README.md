
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Regression with Skew-Normally Distributed Error Term

<!-- badges: start -->

[![](http://cranlogs.r-pkg.org/badges/grand-total/npsf?color=blue)](https://cran.r-project.org/package=npsf)
[![](http://cranlogs.r-pkg.org/badges/last-month/npsf?color=yellow)](https://cran.r-project.org/package=npsf)
[![](https://www.r-pkg.org/badges/version/npsf?color=green)](https://cran.r-project.org/package=npsf)
[![](https://img.shields.io/badge/devel%20version-1.0-green.svg)](https://github.com/OlegBadunenko/snreg)
[![CRAN
checks](https://badges.cranchecks.info/summary/npsf.svg)](https://cran.r-project.org/web/checks/check_results_npsf.html)
[![](https://img.shields.io/github/last-commit/OlegBadunenko/snreg.svg)](https://github.com/OlegBadunenko/snreg/commits/main)

<!-- badges: end -->

<!-- README.md is generated from README.Rmd. Please edit that file -->

The `snreg` package offers a set of methods for conducting regression
analysis when the model errors follow a skew‑normal distribution.

## The framework

The `snreg` package implements the framework developed in

Oleg Badunenko and Daniel J. Henderson (2023). “Production analysis with
asymmetric noise”. Journal of Productivity Analysis, 61(1), 1–18. [DOI
<img src="man/figures/doi.png"  width="12" height="12">](https://doi.org/10.1007/s11123-023-00680-5)

R commands `snreg` and `snsf` estimate models with skew-normal errors
written and maintained by Oleg Badunenko
(<oleg.badunenko@proton.me>).

# 📦 Installing `snreg` R Package from GitHub

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

## Illustration and Uses

This [article](https://olegbadunenko.github.io/snreg/illustration.html)
guides through the code and illustrates the functionality of the package
using

> a subset of the banking data (`banks07`) available in the package.

