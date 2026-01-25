
#' Owen's T Function Variant via C Backend
#'
#' @title Compute Owen-like Function \code{TOwen1(h, a)}
#'
#' @description
#' \code{TOwen1} computes an Owen's \eqn{T}-function variant (or a related
#' special function) for vectors \code{h} and \code{a} by delegating to a compiled
#' C routine named \code{"TOwen1"}. Non-finite inputs in \code{h} or \code{a}
#' yield \code{NA} at the corresponding positions.
#'
#' @details
#' This is a thin R wrapper around a native routine with signature:
#' \preformatted{
#'   void TOwen1(int *n, double *h, double *a, double *out, int *threads)
#' }
#' The C function is expected to fill \code{out[0..(*n-1)]} with results for
#' each pair \code{(h[i], a[i])}. The wrapper:
#' \itemize{
#'   \item Recycles \code{a} to match \code{length(h)} (standard R behavior).
#'   \item Computes values only for finite pairs; non-finite entries become \code{NA}.
#'   \item Forwards a \code{threads} hint to the C implementation.
#' }
#'
#' @section Linking to the C routine:
#' Ensure the symbol \code{"TOwen1"} is available from your package's shared
#' object. With Roxygen2, include (once in your package):
#' \preformatted{
#' #' @useDynLib yourpkg, .registration = TRUE
#' NULL
#' }
#' and register the routine in \code{src/init.c}. If you are not using explicit
#' registration, \code{@useDynLib yourpkg} (without \code{.registration = TRUE})
#' will still allow \code{.C("TOwen1", ...)} provided the symbol is exported.
#'
#' @param h
#' numeric vector of \eqn{h} arguments.
#'
#' @param a
#' numeric vector of \eqn{a} arguments. Must be either the same length as \code{h}
#' or of length 1 (it will be recycled to match \code{h}).
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
#'   # Requires compiled C routine "TOwen1":
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