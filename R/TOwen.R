TOwen <- function(h,a,threads = 1){
  n <- length(h)
  # cat.print(h)
  # cat.print(a)
  # cat.print(threads)
  toInclude <- is.finite(h) & is.finite(a) 
  # toInclude <- rep(TRUE,n)
  if(sum(toInclude) < n){
    tymch1 <- rep(NA, n)
    n1 <- sum(toInclude)
    h1 <- h[toInclude]
    a1 <- a[toInclude]
    tymch1[toInclude] <- .C("TOwen", as.integer(n1), as.double(h1), as.double(a1), TOwenValues = double(n1), as.integer(threads)) $ TOwenValues
  } else {
    tymch1 <- .C("TOwen", as.integer(n), as.double(h), as.double(a), TOwenValues = double(n), as.integer(threads)) $ TOwenValues
  }
  class(tymch1) <- "snreg"
  return(tymch1)
}