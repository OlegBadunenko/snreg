
#' U.S. Commercial Banks Data (2007)
#'
#' @name banks07
#' @docType data
#'
#' @title U.S. Commercial Banks Data
#'
#' @description
#' \code{banks07} is a data frame containing selected variables for
#' 500 U.S. commercial banks, randomly sampled from approximately 5000 banks,
#' based on the dataset of Koetter et al. (2012) for year 2007.
#' The dataset is provided solely for illustration and pedagogical purposes
#' and is not suitable for empirical research.
#'
#' @usage data(banks07)
#'
#' @format
#' A data frame with the following variables:
#' \describe{
#'   \item{\code{year}}{Year (2007).}
#'   \item{\code{id}}{Entity (bank) identifier.}
#'   \item{\code{TA}}{Gross total assets.}
#'   \item{\code{LLP}}{Loan loss provisions.}
#'   \item{\code{Y1}}{Total securities (thousands of USD).}
#'   \item{\code{Y2}}{Total loans and leases (thousands of USD).}
#'   \item{\code{W1}}{Cost of fixed assets divided by the cost of borrowed funds.}
#'   \item{\code{W2}}{Cost of labor (thousands of USD) divided by the cost of borrowed funds.}
#'   \item{\code{ER}}{Equity-to-assets ratio (gross).}
#'   \item{\code{TC}}{Total operating cost.}
#'   \item{\code{LA}}{Ratio of total loans and leases to gross total assets.}
#' }
#'
#' @details
#' The dataset was created by sampling and transforming variables as shown in the
#' section \strong{Examples}.  
#' It is intended to illustrate the usage of functions from this package
#' (e.g. stochastic frontier models with skew-normal noise).
#'
#' @examples
#' \dontrun{
#'
#' data(banks07)
#' head(banks07)
#'
#' dat <- banks07
#' dat$mysample <- TRUE
#'
#' # Example: simple translog cost function
#' spe.tl <-
#'   log(TC) ~ (log(Y1) + log(Y2) + log(W1) + log(W2))^2 +
#'   I(0.5 * log(Y1)^2) + I(0.5 * log(Y2)^2) +
#'   I(0.5 * log(W1)^2) + I(0.5 * log(W2)^2)
#'
#' # Fit a normal-noise MLE model via snreg::lm.mle
#' m.N0 <- snreg::lm.mle(
#'   formula = spe.tl, data = dat, subset = mysample,
#'   ln.var.v = NULL, technique = c("bfgs")
#' )
#'
#' m.N0$value
#' m.N0$coef
#' sqrt(diag(m.N0$vcov))
#'
#' # Fit a skew-normal noise model
#' formSK <- NULL
#' formSV <- NULL
#'
#' init.sk.grid <- seq(0, 3, 0.07)
#' init.sk.grid[1] <- 0.01
#' init.sk.grid <- sort(c(-init.sk.grid, init.sk.grid))
#'
#' trial <- NULL
#' for (in.sk in init.sk.grid) {
#'   m.temp <- snreg::snreg(
#'     formula = spe.tl, data = dat, subset = mysample,
#'     start.sk = in.sk, ln.var.v = formSV, skew.v = formSK,
#'     technique = c("bfgs")
#'   )
#'   trial <- rbind(trial, c(in.sk, m.temp$value))
#' }
#'
#' best.sk <- trial[order(trial[, 2]),][nrow(trial), 1]
#'
#' m.SN0 <- snreg::snreg(
#'   formula = spe.tl, data = dat, subset = mysample,
#'   start.sk = best.sk, ln.var.v = formSV, skew.v = formSK,
#'   technique = c("nr")
#' )
#'
#' # Fit a stochastic frontier model with exponential inefficiency
#' m.SF1 <- snreg::snsf(
#'   formula = spe.tl, data = dat, subset = mysample,
#'   distribution = "e", prod = FALSE,
#'   start.sk = -1, ln.var.u = NULL, ln.var.v = NULL, skew.v = NULL,
#'   technique = c("bfgs")
#' )
#'
#' }
#'
#' @source
#' \url{http://qed.econ.queensu.ca/jae/2014-v29.2/restrepo-tobon-kumbhakar/}
#'
#' @references
#' Koetter, M., Kolari, J., & Spierdijk, L. (2012).
#' \emph{Enjoying the quiet life under deregulation? Evidence from adjusted Lerner indices for U.S. banks}.
#' Review of Economics and Statistics, \bold{94}(2), 462–480.
#'
#' Restrepo-Tobon, D. & Kumbhakar, S. (2014).
#' \emph{Enjoying the quiet life under deregulation? Not Quite}.
#' Journal of Applied Econometrics, \bold{29}(2), 333–343.
#'
#' @keywords datasets
NULL
