coef.snreg <- function( obj, ... ) {
  if(is.null(obj$coef) | is.null(obj)){
    stop( paste("No results are available in ",obj,"") )
  }
  return( obj$coef )
}