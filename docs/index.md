# Regression with Skew-Normally Distributed Error Term

The `snreg` package offers a set of methods for conducting regression
analysis when the model errors follow a skew‑normal distribution.

## The framework

The `snreg` package implements the framework developed in

Oleg Badunenko and Daniel J. Henderson (2023). “Production analysis with
asymmetric noise”. Journal of Productivity Analysis, 61(1), 1–18. [DOI
![](reference/figures/doi.png)](https://doi.org/10.1007/s11123-023-00680-5)

R commands `snreg` and `snsf` estimate models with skew-normal errors
written and maintained by Oleg Badunenko
(<oleg.badunenko@brunel.ac.uk>).

## Acknowledgments

The R package `snreg` computes Owen’s *T* function using C code written
by John Burkardt. This implementation, distributed under the MIT
license, is publicly accessible at
<https://people.sc.fsu.edu/~jburkardt/c_src/owen/owen.html>.

# 📦 Installing `snreg` R Package from GitHub

### 1. Install the `devtools` package

> Install `devtools` from CRAN (if you haven’t already):

``` r

# Install only if not already installed
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
```

### 2. Install the package from GitHub

> Install the `snreg` package from GitHub. In this code, you are
> installing the `snreg` package created by `OlegBadunenko`.

``` r

# Install only if not already installed
if (!requireNamespace("snreg", quietly = TRUE)) {
  remotes::install_github("OlegBadunenko/snreg", dependencies = TRUE, build_vignettes = FALSE)
} else {
  message("Package 'snreg' is already installed; skipping.")
}
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

This
[article](https://olegbadunenko.github.io/snreg/articles/illustration.html)
guides through the code and illustrates the functionality of the package
using

> a subset of the banking data (`banks07`) available in the package.
