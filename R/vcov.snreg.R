vcov.snreg <- function( obj, ... ) {
  if(is.null(obj$vcov) | is.null(obj)){
    stop( paste("No vcov are available in ",obj,"") )
  }
  return( obj$vcov )
}