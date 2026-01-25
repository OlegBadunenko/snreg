
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
#' ## ------------------------------------------------------------------
#' ## Construct sample panel dataset (banks00_07)
#' ## ------------------------------------------------------------------
#'
#' # Download data from the link in "Source"
#' banks00_07 <- read.delim("2b_QLH.txt")
#'
#' # rename 'entity' to 'id'
#' colnames(banks00_07)[colnames(banks00_07) == "entity"] <- "id"
#'
#' # keep only years 2000–2007
#' banks00_07 <- banks00_07[
#'   banks00_07$year >= 2000 & banks00_07$year <= 2007, ]
#'
#' # restrict sample to interquartile range of total assets
#' q1q3 <- quantile(banks00_07$TA, probs = c(.25, .75))
#' banks00_07 <- banks00_07[
#'   banks00_07$TA >= q1q3[1] & banks00_07$TA <= q1q3[2], ]
#'
#' # generate required variables
#' banks00_07$TC <- banks00_07$TOC
#' banks00_07$ER <- banks00_07$Z  / banks00_07$TA   # Equity ratio
#' banks00_07$LA <- banks00_07$Y2 / banks00_07$TA   # Loans-to-assets ratio
#'
#' # keep only needed variables
#' keep.vars <- c("id", "year", "Ti", "TC", "Y1", "Y2", "W1","W2",
#'                "ER", "LA", "TA", "LLP")
#' banks00_07 <- banks00_07[, colnames(banks00_07) %in% keep.vars]
#'
#' # number of periods per id
#' t0 <- as.vector( by(banks00_07$id, banks00_07$id,
#'                     FUN = function(qq) length(qq)) )
#' banks00_07$Ti <- rep(t0, times = t0)
#'
#' # keep if Ti > 4
#' banks00_07 <- banks00_07[banks00_07$Ti > 4, ]
#'
#' # complete observations only
#' banks00_07 <- banks00_07[complete.cases(banks00_07), ]
#'
#' # sample 500 banks at random
#' set.seed(816376586)
#' id_names <- unique(banks00_07$id)
#' ids2choose <- sample(id_names, 500)
#' banks00_07 <- banks00_07[banks00_07$id %in% ids2choose, ]
#'
#' # recompute Ti
#' t0 <- as.vector( by(banks00_07$id, banks00_07$id,
#'                     FUN = function(qq) length(qq)) )
#' banks00_07$Ti <- rep(t0, times = t0)
#' banks00_07 <- banks00_07[banks00_07$Ti > 4, ]
#'
#' # sort
#' banks00_07 <- banks00_07[order(banks00_07$id, banks00_07$year), ]
#' 
#' 
# keep only year 2007
#' banks07 <- banks00_07[banks00_07$year == 2007, ]
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
