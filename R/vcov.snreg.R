
#' Variance–Covariance Matrix for snreg Objects
#'
#' @title Extract the Variance–Covariance Matrix
#'
#' @description
#' \code{vcov.snreg} is the \code{vcov} S3 method for objects of class \code{"snreg"}.
#' It returns the model-based variance–covariance matrix stored in the fitted object.
#'
#' @param obj
#' an object of class \code{"snreg"}, typically returned by \code{\link{snreg}}.
#'
#' @param ...
#' additional arguments (currently unused).
#'
#' @return
#' A numeric matrix containing the variance–covariance of the estimated parameters.
#'
#' @details
#' This method expects a fitted \code{"snreg"} object.
#' 
#' This method simply returns the \code{vcov} component stored in \code{obj}.
#' If your estimator did not compute standard errors (e.g., because estimation
#' hasn’t been run yet in a scaffold), this field may be \code{NULL}, and the
#' method will error accordingly.
#'
#' @seealso
#' \code{\link{snreg}}, \code{\link{summary.snreg}}
#'
#' @examples
#' \dontrun{
#'   library(snsf)
#'
#'   data("banks07")
#'   head(banks07)
#'   # V <- vcov(m)
#'   # diag(V)
#' }
#'
#' @export
vcov.snreg <- function(obj, ...) {
  if (is.null(obj)) {
    stop("Argument 'obj' is NULL; expected a fitted 'snreg' object.", call. = FALSE)
  }
  if (is.null(obj$vcov)) {
    stop("No variance–covariance matrix available in 'obj$vcov'.", call. = FALSE)
  }
  obj$vcov
}


# vcov.snreg <- function( obj, ... ) {
#   if(is.null(obj$vcov) | is.null(obj)){
#     stop( paste("No vcov are available in ",obj,"") )
#   }
#   return( obj$vcov )
# }