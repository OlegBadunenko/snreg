su <- function(x, mat.var.in.col = TRUE, digits = 4, probs = c(0.1, 0.25, 0.5, 0.75, 0.9), print = FALSE){

  xvec2 <- xvec1 <- FALSE

  if(is.matrix(x)){
    if(min(dim(x)) == 1){
      # xvec1 <- TRUE
      if(which(dim(x) == 1) == 2){
        mynames <- colnames(x)
      } else {
        mynames <- rownames(x)
        x <- t(x)
      }
      # x <- as.vector(x)
    } else {
      if(!mat.var.in.col){
        x <- t(x)
      }
      mynames <- colnames(x)
    }
    # print(mynames)
    if(is.null(mynames)) mynames <- paste("Var", seq_len(ncol(x)), sep = "")
    x <- as.data.frame(x)
    # print(x)
    # mynames <- colnames(x)
  } # end if matrix

  if(is.vector(x)){
    xvec2 <- TRUE
    mynames <- deparse(substitute(x))
    x <- data.frame(Var1 = x)
  } # end if vector

  # cat("nymanes", sep ="")
  # print(mynames)

  if(!is.vector(x) & !is.matrix(x) & !is.data.frame(x)){
    stop("Provide vector, matrix, or data.frame")
  } else {
    t1 <- apply(x, 2, function(x) c(Obs = length(x), NAs = sum(is.na(x)), Mean = mean(x, na.rm = TRUE), StDev = sd(x, na.rm = TRUE), IQR = IQR(x, na.rm = TRUE), Min = min(x, na.rm = TRUE), quantile( x, probs = probs, na.rm = TRUE ), Max = max(x, na.rm = TRUE)))
    # print(t1)
    # print(mynames)
    # print(class(t1))
    # print(dim(t1))
    if(xvec2 & !xvec1) colnames(t1) <- mynames
    if(print) print(t(t1), digits = digits)
  }
  tymch <- t(t1)
  class(tymch) <- "snreg"
  return(tymch)
}