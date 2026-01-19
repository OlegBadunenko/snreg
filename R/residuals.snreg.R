residuals.snreg <- function( obj, esample = TRUE, ... ) {
  if(is.null(obj$resid) | is.null(obj)){
    stop( paste("No residuals are available in ",obj,"") )
  }
  if(esample){
    my.res <- rep(NA, length(obj$esample))
    my.res[obj$esample] <- obj$resid
  } else {
    my.res <- obj$resid
  }
  return( my.res )
}