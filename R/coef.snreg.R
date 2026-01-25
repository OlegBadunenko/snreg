
#' Coefficients from an snreg Model
#'
#' @title Extract Model Coefficients
#'
#' @description
#' \code{coef.snreg} is the S3 method for extracting the estimated regression
#' coefficients from an object of class \code{"snreg"}.
#'
#' @param obj
#' an object of class \code{"snreg"}, typically returned by \code{\link{snreg}}.
#'
#' @param ...
#' additional arguments (currently unused).
#'
#' @return
#' A numeric vector containing the model coefficients.
#'
#' @details
#' This method simply returns the \code{coef} component stored inside the fitted
#' \code{"snreg"} object. If the object does not contain coefficient estimates
#' (e.g., if estimation was not completed in a scaffold), an informative error
#' is raised.
#'
#' @examples
#' \dontrun{
#'   m <- snreg(y ~ x1 + x2, data = df)
#'   coef(m)
#' }
#'
#' @export
coef.snreg <- function(obj, ...) {
  if (is.null(obj)) {
    stop("Argument 'obj' is NULL; expected a fitted 'snreg' object.", call. = FALSE)
  }
  if (is.null(obj$coef)) {
    stop("No coefficients available in 'obj$coef'.", call. = FALSE)
  }
  obj$coef
}

# coef.snreg <- function( obj, ... ) {
#   if(is.null(obj$coef) | is.null(obj)){
#     stop( paste("No results are available in ",obj,"") )
#   }
#   return( obj$coef )
# }