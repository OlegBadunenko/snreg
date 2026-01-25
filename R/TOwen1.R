
#' Owen's T Function Variant via C Backend
#'
#' @title Compute Owen-like Function \code{TOwen1(h, a)}
#'
#' @description
#' \code{TOwen1} computes an Owen's \eqn{T}-function variant (or a related
#' special function) for vectors \code{h} and \code{a} based on the \code{tha} function in \url{https://people.sc.fsu.edu/~jburkardt/c_src/owen/owen.html}. Non-finite inputs in \code{h} or \code{a}
#' yield \code{NA} at the corresponding positions.
#'
#' @details
#' This is a thin R wrapper around a native routine with signature:
#' \preformatted{
#'   void TOwen1(int *n, double *h, double *a, double *out, int *threads)
#' }
#'
#' @param threads
#' integer. Number of threads to request from the C implementation (if supported).
#' Default is \code{1}.
#'
#' @return
#' A numeric vector of length \code{length(h)} with the computed values. Elements
#' where either \code{h} or \code{a} is non-finite are \code{NA}. The returned
#' vector is given class \code{"snreg"} for downstream compatibility.
#'
#' @examples
#' \dontrun{
#'   library(snreg)
#' 
#'   # Basic usage. Vectorized 'a':
#'   h <- c(-1, 0, 1, 2)
#'   a <- 0.3
#'   TOwen1(h, a)
#'
#'   # Vectorized 'a' with non-finite entries:
#'   a2 <- c(0.2, NA, 1, Inf)
#'   TOwen1(h, a2)
#' }
#'
#' @seealso
#' \code{\link{TOwen}}
#'
#' @export
TOwen1 <- function(h, a, threads = 1) {
  n <- length(h)
  
  # Recycle 'a' if needed to match 'h'
  if (length(a) != n) {
    a <- rep(a, length.out = n)
  }
  
  # Compute only for finite pairs (h, a); others become NA
  toInclude <- is.finite(h) & is.finite(a)
  
  if (sum(toInclude) < n) {
    out <- rep(NA_real_, n)
    n1  <- sum(toInclude)
    if (n1 > 0) {
      h1 <- as.double(h[toInclude])
      a1 <- as.double(a[toInclude])
      out[toInclude] <- .C(
        "TOwen1",
        as.integer(n1),
        h1,
        a1,
        TOwenValues = double(n1),
        as.integer(threads)
      )$TOwenValues
    }
  } else {
    out <- .C(
      "TOwen1",
      as.integer(n),
      as.double(h),
      as.double(a),
      TOwenValues = double(n),
      as.integer(threads)
    )$TOwenValues
  }
  
  class(out) <- "snreg"
  out
}


# OLD before AI

# TOwen1 <- function(h,a,threads = 1){
#   n <- length(h)
#   # cat.print(h)
#   # cat.print(a)
#   # cat.print(threads)
#   toInclude <- is.finite(h) & is.finite(a) 
#   # toInclude <- rep(TRUE,n)
#   if(sum(toInclude) < n){
#     tymch1 <- rep(NA, n)
#     n1 <- sum(toInclude)
#     h1 <- h[toInclude]
#     a1 <- a[toInclude]
#     tymch1[toInclude] <- .C("TOwen1", as.integer(n1), as.double(h1), as.double(a1), TOwenValues = double(n1), as.integer(threads)) $ TOwenValues
#   } else {
#     tymch1 <- .C("TOwen1", as.integer(n), as.double(h), as.double(a), TOwenValues = double(n), as.integer(threads)) $ TOwenValues
#   }
#   class(tymch1) <- "snreg"
#   return(tymch1)
# }