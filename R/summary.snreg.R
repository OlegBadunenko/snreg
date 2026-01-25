
#' Summary Method for snreg Objects
#'
#' @title Summary for Skew-Normal Regression Models
#'
#' @description
#' Produces a summary object for objects of class \code{"snreg"}.
#' The function assigns the class \code{"summary.snreg"} to the fitted model
#' object, enabling a dedicated print method (\code{print.summary.snreg}) to
#' display results in a structured format.
#'
#' @param obj
#' an object of class \code{"snreg"}, typically returned by \code{\link{snreg}}.
#'
#' @param ...
#' additional arguments (currently not used).
#'
#' @details
#' \code{summary.snreg} does not modify the contents of the object; it only
#' updates the class attribute to \code{"summary.snreg"}. The corresponding
#' print method (\code{\link{print.summary.snreg}}) is responsible for
#' formatting and displaying estimation details, such as convergence criteria,
#' log-likelihood, coefficient tables, and (if present) heteroskedastic and
#' skewness components. The print method shown below assumes the presence of
#' some internal helper functions:
#' \itemize{
#'   \item \code{.timing(time_object, prefix)} — to print elapsed time,
#'   \item \code{.printgtresfhet(res, digits, Kb, Kv0, Ku0, Kvi, Kui, na.print, max.name.length)} — to format coefficient blocks,
#'   \item \code{.su(x, ...)} — to summarize vectors/matrices, used for efficiency summaries.
#' }
#' Ensure these are defined in your package if you wish to use this print method unchanged.
#'
#' @return
#' An object of class \code{"summary.snreg"}, identical to the input \code{obj}
#' except for its class attribute.
#'
#' @seealso
#' \code{\link{snreg}}, \code{\link{print.summary.snreg}}
#'
#' @examples
#' \dontrun{
#'   # Example (assuming 'banks07' is available and snreg is fitted accordingly)
#'   # m <- snreg(TC ~ Y1 + Y2, data = banks07)
#'   # s <- summary(m)
#'   # print(s)
#' }
#'
#' @export
summary.snreg <- function( obj, ...) {
  class(obj) <- "summary.snreg"
  return(obj)
}


#' Print Method for Summary of snreg Objects
#'
#' @title Print Summary of snreg Results
#'
#' @description
#' Prints the contents of a \code{"summary.snreg"} object in a structured
#' format. The method reports convergence status (based on gradient-Hessian
#' scaling), log-likelihood, estimation results, and—when present—summaries
#' for technical/cost efficiencies and marginal effects.
#'
#' @param obj
#' an object of class \code{"summary.snreg"} (produced by \code{\link{summary.snreg}}).
#'
#' @param digits
#' integer indicating the number of digits to print; default \code{NULL}
#' (internally set to 4).
#'
#' @param ...
#' additional arguments (currently unused).
#'
#' @details
#' This method expects the input object \code{obj} to contain fields produced by
#' your estimation routine, including (but not limited to):
#' \code{results}, \code{gHg}, \code{lmtol}, \code{LM}, \code{bhhh}, \code{gg},
#' \code{esttime}, \code{ll}, \code{prod}, \code{toinclude}, \code{Kv0}, \code{Ku0},
#' \code{Kvi}, \code{Kui}, \code{Kb}, \code{nt}, \code{n}, \code{dat.descr},
#' and efficiency objects such as \code{te_i0}, \code{te_it}, \code{te_over},
#' \code{ce_i0}, \code{ce_it}, \code{ce_over}, as well as marginal effects summaries
#' \code{me_i0}, \code{me_it} and their elasticities.
#'
#' It also uses internal helpers \code{.timing}, \code{.printgtresfhet},
#' and \code{.su}. Make sure those are available in your package namespace.
#'
#' @return
#' The input \code{obj} is returned (invisibly) after printing.
#'
#' @seealso
#' \code{\link{summary.snreg}}
#'
#' @export
print.summary.snreg <- function( obj, digits = NULL, ... ) {
  # Defensive: handle empty results gracefully when computing max.name.length
  rn <- tryCatch(row.names(obj$results), error = function(e) NULL)
  max.name.length <- if (!is.null(rn) && length(rn)) max(nchar(rn)) else 12
  
  if( is.null( digits ) ) {
    digits <- 4
  }
  if(obj$gHg < obj$lmtol){
    cat("\nConvergence given g*inv(H)*g' = ",
        formatC(obj$gHg, digits = 1, format = "e"),
        " < lmtol(", obj$lmtol, ")\n", sep = "")
    # cat(" was reached in ",obj$counts[2]," iteration(s)\n", sep = "")
  }
  else {
    cat("\nCriterion g*inv(H)*g' = ",
        formatC(obj$LM, digits = 1, format = "e"),
        " > lmtol(", obj$lmtol, ")\n", sep = "")
    if(obj$bhhh){
      cat("Note that Hessian is computed as outer product (BHHH)\n", sep = "")
      cat("Criterion g'g = ", obj$gg, "\n", sep = "")
    }
    warning("Convergence given g*inv(H)*g' is still not reached; one of optim's convergence criteria is used",
            call. = FALSE)
  }
  .timing(obj$esttime, "Log likelihood maximization completed in ")
  cat("Log likelihood = ", formatC(obj$ll, digits = 4, format = "f"), "\n",  sep = "")
  # cat("____________________________________________________\n")
  cat("", rep("_", max.name.length + 42 - 1), "", "\n", sep = "")
  
  cat("\nStochastic ", ifelse(obj$prod, "production", "cost"), " frontier   model\n",  sep = "")
  cat("\nDistributional assumptions\n\n", sep = "")
  Assumptions <- rep("heteroskedastic", 4)
  if(obj$Kv0 == 1) Assumptions[1] <- "homoskedastic"
  if(obj$Ku0 == 1) Assumptions[2] <- "homoskedastic"
  if(obj$Kvi == 1) Assumptions[3] <- "homoskedastic"
  if(obj$Kui == 1) Assumptions[4] <- "homoskedastic"
  a1 <- data.frame(
    Component = c("Random effects:", "Persistent ineff.: ", "Random noise:", "Transient ineff.: "),
    Distribution = c("normal", "half-normal  ", "normal", "half-normal  "),
    Assumption = Assumptions
  )
  print(a1[obj$toinclude, ], quote = FALSE, right = FALSE)
  
  est.spd.left  <- floor( (max.name.length + 42 - 29) / 2 )
  est.spd.right <- max.name.length + 42 - 29 - est.spd.left
  
  cat("\n", rep("-", est.spd.left), " Summary of the panel data: ", rep("-", est.spd.right), "\n\n", sep ="")
  cat("   Number of obs       (NT) =", obj$nt, "", "\n")
  cat("   Number of groups     (N) =", obj$n, "", "\n")
  cat("   Obs per group: (T_i) min =", obj$dat.descr[3], "", "\n")
  cat("                        avg =", obj$dat.descr[4], "", "\n")
  cat("                        max =", obj$dat.descr[5], "", "\n")
  
  est.rez.left  <- floor( (max.name.length + 42 - 22) / 2 )
  est.rez.right <- max.name.length + 42 - 22 - est.rez.left
  cat("\n", rep("-", est.rez.left), " Estimation results: ", rep("-", est.rez.right), "\n\n", sep ="")
  
  .printgtresfhet(obj$results, digits = digits, obj$Kb, obj$Kv0, obj$Ku0, obj$Kvi, obj$Kui,
                  na.print = "NA", max.name.length)
  
  if(obj$prod){
    myeff <- "technical"
    ndots <- "...."
    e_i0  <- obj$te_i0
    e_it  <- obj$te_it
    e_ov  <- obj$te_over
  }
  else {
    myeff <- "cost"
    ndots <- "........."
    e_i0  <- obj$ce_i0
    e_it  <- obj$ce_it
    e_ov  <- obj$ce_over
  }
  
  sum.te <- paste0(" Summary of ", myeff, " efficiencies: ")
  est.cle.left  <- floor( (max.name.length + 42 - nchar(sum.te) - 1) / 2 )
  est.cle.right <- max.name.length + 42 - nchar(sum.te) - 1 - est.cle.left
  if(est.cle.left  <= 0) est.cle.left  <- 1
  if(est.cle.right <= 0) est.cle.right <- 1
  cat("\n", rep("-", est.cle.left), sum.te, rep("-", est.cle.right), "\n\n", sep ="")
  .su(list(e_i0, e_it, e_ov), print = TRUE, width = 5, format = "fg",
      drop0trailing = FALSE, names = c("Persistent", "Transient", "Overall"))
  
  # ME
  if(obj$Ku0 > 1){
    cat("", rep("_", max.name.length + 42 - 1), "", "\n", sep = "")
    
    sum.me_i0 <- paste0(" Summary of ME of Zi0 on persistent ", myeff, " eff-s: ")
    est.cle.left  <- floor( (max.name.length + 42 - nchar(sum.me_i0) - 1) / 2 )
    est.cle.right <- max.name.length + 42 - nchar(sum.me_i0) - 1 - est.cle.left
    if(est.cle.left  <= 0) est.cle.left  <- 2
    if(est.cle.right <= 0) est.cle.right <- 2
    cat("\n", rep("-", est.cle.left), sum.me_i0, rep("-", est.cle.right), "\n\n", sep ="")
    .su(obj$me_i0, print = TRUE)
    
    sum.me_i0 <- paste0(" Summary of elasticity of E(u_i0) w.r.t Zi0: ")
    est.cle.left  <- floor( (max.name.length + 42 - nchar(sum.me_i0) - 1) / 2 )
    est.cle.right <- max.name.length + 42 - nchar(sum.me_i0) - 1 - est.cle.left
    if(est.cle.left  <= 0) est.cle.left  <- 2
    if(est.cle.right <= 0) est.cle.right <- 2
    cat("\n", rep("-", est.cle.left), sum.me_i0, rep("-", est.cle.right), "\n\n", sep ="")
    .su(obj$me_i0_elast, print = TRUE)
  }
  if(obj$Kui > 1){
    cat("", rep("_", max.name.length + 42 - 1), "", "\n", sep = "")
    
    sum.me_i0 <- paste0(" Summary of ME of Zit on transient ", myeff, " eff-s: ")
    est.cle.left  <- floor( (max.name.length + 42 - nchar(sum.me_i0) - 1) / 2 )
    est.cle.right <- max.name.length + 42 - nchar(sum.me_i0) - 1 - est.cle.left
    if(est.cle.left  <= 0) est.cle.left  <- 2
    if(est.cle.right <= 0) est.cle.right <- 2
    cat("\n", rep("-", est.cle.left), sum.me_i0, rep("-", est.cle.right), "\n\n", sep ="")
    .su(obj$me_it, print = TRUE)
    
    sum.me_i0 <- paste0(" Summary of elasticity of E(u_it) w.r.t Zit: ")
    est.cle.left  <- floor( (max.name.length + 42 - nchar(sum.me_i0) - 1) / 2 )
    est.cle.right <- max.name.length + 42 - nchar(sum.me_i0) - 1 - est.cle.left
    if(est.cle.left  <= 0) est.cle.left  <- 2
    if(est.cle.right <= 0) est.cle.right <- 2
    cat("\n", rep("-", est.cle.left), sum.me_i0, rep("-", est.cle.right), "\n\n", sep ="")
    .su(obj$me_it_elast, print = TRUE)
  }
  
  invisible( obj )
}
