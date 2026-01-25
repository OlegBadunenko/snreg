
#' Compact Summary Statistics for Vectors, Matrices, and Data Frames
#'
#' @title Summary Utility (\code{su})
#'
#' @description
#' Computes a compact table of summary statistics for each variable in a vector,
#' matrix, or data frame. The following metrics are returned per variable:
#' number of observations (\code{Obs}), missing values (\code{NAs}), mean,
#' standard deviation (\code{StDev}), interquartile range (\code{IQR}),
#' minimum (\code{Min}), user-specified quantiles (\code{probs}), and maximum (\code{Max}).
#'
#' @param x
#' a numeric vector, matrix, or data frame. For matrices, variables are assumed
#' to be in columns; set \code{mat.var.in.col = FALSE} to treat rows as variables.
#'
#' @param mat.var.in.col
#' logical. If \code{TRUE} (default), a matrix is interpreted as variables in columns.
#' If \code{FALSE}, the matrix is transposed so that rows are treated as variables.
#'
#' @param digits
#' integer. Number of digits to use when printing (only affects printed output when
#' \code{print = TRUE}). Default is \code{4}.
#'
#' @param probs
#' numeric vector of probabilities in \eqn{[0, 1]} for which quantiles are computed.
#' Default is \code{c(0.1, 0.25, 0.5, 0.75, 0.9)}.
#'
#' @param print
#' logical. If \code{TRUE}, prints the transposed summary table using the specified
#' number of digits. Default is \code{FALSE}.
#'
#' @details
#' Input handling:
#' \itemize{
#'   \item If \code{x} is a matrix with a single row or column, it is treated like a vector.
#'         Column or row names are used (if available). Otherwise, a default name is created.
#'   \item If \code{x} is a matrix with multiple variables, variables are taken as columns.
#'         Use \code{mat.var.in.col = FALSE} to transpose and treat rows as variables.
#'   \item If \code{x} is a vector, its deparsed symbol name is used as the variable name.
#'   \item If \code{x} is a data frame, each column is summarized.
#' }
#'
#' Missing values are excluded in all summary computations.
#'
#' @return
#' A matrix (coercible to \code{data.frame}) where each row corresponds to a variable
#' and columns contain the summary statistics:
#' \code{Obs}, \code{NAs}, \code{Mean}, \code{StDev}, \code{IQR}, \code{Min},
#' the requested \code{probs} quantiles (named), and \code{Max}.  
#' The returned object is given class \code{"snreg"} for compatibility with
#' package-specific print/summarization methods.
#'
#' @examples
#' \dontrun{
#'   # Vector
#'   set.seed(1)
#'   v <- rnorm(100)
#'   su(v, print = TRUE)
#'
#'   # Matrix: variables in columns
#'   M <- cbind(x = rnorm(50), y = runif(50))
#'   su(M)
#'
#'   # Matrix: variables in rows
#'   Mr <- rbind(x = rnorm(50), y = runif(50))
#'   su(Mr, mat.var.in.col = FALSE)
#'
#'   # Data frame
#'   DF <- data.frame(a = rnorm(30), b = rexp(30), c = rbinom(30, 1, 0.3))
#'   out <- su(DF)
#'   head(out)
#' }
#'
#' @export
su <- function(x, mat.var.in.col = TRUE, digits = 4,
               probs = c(0.1, 0.25, 0.5, 0.75, 0.9), print = FALSE) {
  
  xvec2 <- xvec1 <- FALSE
  
  if (is.matrix(x)) {
    if (min(dim(x)) == 1) {
      # effectively a vector held as a matrix
      if (which(dim(x) == 1) == 2) {
        mynames <- colnames(x)
      } else {
        mynames <- rownames(x)
        x <- t(x)
      }
    } else {
      if (!isTRUE(mat.var.in.col)) {
        x <- t(x)
      }
      mynames <- colnames(x)
    }
    if (is.null(mynames)) mynames <- paste("Var", seq_len(ncol(x)), sep = "")
    x <- as.data.frame(x)
  }
  
  if (is.vector(x)) {
    xvec2 <- TRUE
    mynames <- deparse(substitute(x))
    x <- data.frame(Var1 = x)
  }
  
  if (!is.vector(x) && !is.matrix(x) && !is.data.frame(x)) {
    stop("Provide vector, matrix, or data.frame", call. = FALSE)
  } else {
    t1 <- apply(
      x, 2,
      function(x) c(
        Obs   = length(x),
        NAs   = sum(is.na(x)),
        Mean  = mean(x, na.rm = TRUE),
        StDev = stats::sd(x, na.rm = TRUE),
        IQR   = stats::IQR(x, na.rm = TRUE),
        Min   = min(x, na.rm = TRUE),
        stats::quantile(x, probs = probs, na.rm = TRUE),
        Max   = max(x, na.rm = TRUE)
      )
    )
    if (xvec2 && !xvec1) colnames(t1) <- mynames
    if (isTRUE(print)) utils::capture.output(print(t(t1), digits = digits), type = "output") |> cat(sep = "\n")
  }
  tymch <- t(t1)
  class(tymch) <- "snreg"
  return(tymch)
}

# OLD before AI
# su <- function(x, mat.var.in.col = TRUE, digits = 4, probs = c(0.1, 0.25, 0.5, 0.75, 0.9), print = FALSE){
# 
#   xvec2 <- xvec1 <- FALSE
# 
#   if(is.matrix(x)){
#     if(min(dim(x)) == 1){
#       # xvec1 <- TRUE
#       if(which(dim(x) == 1) == 2){
#         mynames <- colnames(x)
#       } else {
#         mynames <- rownames(x)
#         x <- t(x)
#       }
#       # x <- as.vector(x)
#     } else {
#       if(!mat.var.in.col){
#         x <- t(x)
#       }
#       mynames <- colnames(x)
#     }
#     # print(mynames)
#     if(is.null(mynames)) mynames <- paste("Var", seq_len(ncol(x)), sep = "")
#     x <- as.data.frame(x)
#     # print(x)
#     # mynames <- colnames(x)
#   } # end if matrix
# 
#   if(is.vector(x)){
#     xvec2 <- TRUE
#     mynames <- deparse(substitute(x))
#     x <- data.frame(Var1 = x)
#   } # end if vector
# 
#   # cat("nymanes", sep ="")
#   # print(mynames)
# 
#   if(!is.vector(x) & !is.matrix(x) & !is.data.frame(x)){
#     stop("Provide vector, matrix, or data.frame")
#   } else {
#     t1 <- apply(x, 2, function(x) c(Obs = length(x), NAs = sum(is.na(x)), Mean = mean(x, na.rm = TRUE), StDev = sd(x, na.rm = TRUE), IQR = IQR(x, na.rm = TRUE), Min = min(x, na.rm = TRUE), quantile( x, probs = probs, na.rm = TRUE ), Max = max(x, na.rm = TRUE)))
#     # print(t1)
#     # print(mynames)
#     # print(class(t1))
#     # print(dim(t1))
#     if(xvec2 & !xvec1) colnames(t1) <- mynames
#     if(print) print(t(t1), digits = digits)
#   }
#   tymch <- t(t1)
#   class(tymch) <- "snreg"
#   return(tymch)
# }