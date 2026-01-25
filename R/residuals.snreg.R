
#' Extract Residuals from an snreg Model
#'
#' @title Residuals for snreg Objects
#'
#' @description
#' \code{residuals.snreg} is the S3 method for extracting residuals from a fitted
#' \code{snreg} model. Residuals may be returned either for the full data or only
#' for the estimation sample.
#'
#' @param obj
#' an object of class \code{"snreg"}, typically produced by \code{\link{snreg}}.
#'
#' @param esample
#' logical. If \code{TRUE} (default), residuals are returned only for observations
#' used in estimation (others are \code{NA}).  
#' If \code{FALSE}, the raw vector of residuals (\code{obj$resid}) is returned.
#'
#' @param ...
#' additional arguments (currently unused).
#'
#' @return
#' A numeric vector of residuals.  
#' If \code{esample = TRUE}, the vector matches the length of the original data
#' and contains \code{NA} for non-estimation observations.  
#' If \code{esample = FALSE}, only the computed residuals are returned.
#'
#' @details
#' This method simply accesses the \code{obj$resid} component of a fitted
#' \code{"snreg"} object. An informative error is produced if residuals are not
#' available.
#'
#' @seealso
#' \code{\link{snreg}}, \code{\link{fitted.snreg}}, \code{\link{coef.snreg}}
#'
#' @examples
#' \dontrun{
#'   m <- snreg(y ~ x1 + x2, data = df)
#'
#'   # Residuals for estimation sample only
#'   residuals(m)
#'
#'   # Residuals for all observations
#'   residuals(m, esample = FALSE)
#' }
#'
#' @export
residuals.snreg <- function(obj, esample = TRUE, ...) {
  
  if (is.null(obj)) {
    stop("Argument 'obj' is NULL; expected a fitted 'snreg' object.", call. = FALSE)
  }
  if (is.null(obj$resid)) {
    stop("No residuals are available in 'obj$resid'.", call. = FALSE)
  }
  
  if (isTRUE(esample)) {
    if (is.null(obj$esample)) {
      stop("Object does not contain an 'esample' indicator.", call. = FALSE)
    }
    out <- rep(NA_real_, length(obj$esample))
    out[obj$esample] <- obj$resid
    return(out)
  }
  
  # esample = FALSE → return raw residuals
  obj$resid
}


# residuals.snreg <- function( obj, esample = TRUE, ... ) {
#   if(is.null(obj$resid) | is.null(obj)){
#     stop( paste("No residuals are available in ",obj,"") )
#   }
#   if(esample){
#     my.res <- rep(NA, length(obj$esample))
#     my.res[obj$esample] <- obj$resid
#   } else {
#     my.res <- obj$resid
#   }
#   return( my.res )
# }