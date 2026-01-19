print.summary.snreg <- function( obj, digits = NULL, ... ) {
  max.name.length <- max(nchar(row.names(obj$results)))
  if( is.null( digits ) ) {
    digits <- 4
  }
  if(obj$gHg < obj$lmtol){
    cat("\nConvergence given g*inv(H)*g' = ",formatC(obj$gHg,digits=1,format="e")," < lmtol(",obj$lmtol,")\n", sep = "")
    # cat(" was reached in ",obj$counts[2]," iteration(s)\n", sep = "")
  }
  else {
    cat("\nCriterion g*inv(H)*g' = ",formatC(obj$LM,digits=1,format="e")," > lmtol(",obj$lmtol,")\n", sep = "")
    if(obj$bhhh){
      cat("Note that Hessian is computed as outer product (BHHH)\n", sep = "")
      cat("Criterion g'g = ",obj$gg,"\n", sep = "")
    }
    warning("Convergence given g*inv(H)*g' is still not reached; one of optim's convergence criteria is used", call. = FALSE)
  }
  .timing(obj$esttime, "Log likelihood maximization completed in ")
  cat("Log likelihood = ",formatC(obj$ll,digits=4,format="f"),"\n",  sep = "")
  # cat("____________________________________________________\n")
  cat("",rep("_", max.name.length+42-1),"", "\n", sep = "")

  cat("\nStochastic ",ifelse(obj$prod, "production", "cost")," frontier   model\n",  sep = "")
  cat("\nDistributional assumptions\n\n", sep = "")
  Assumptions <- rep("heteroskedastic",4)
  if(obj$Kv0==1) Assumptions[1] <- "homoskedastic"
  if(obj$Ku0==1) Assumptions[2] <- "homoskedastic"
  if(obj$Kvi==1) Assumptions[3] <- "homoskedastic"
  if(obj$Kui==1) Assumptions[4] <- "homoskedastic"
  a1 <- data.frame(
    Component = c("Random effects:","Persistent ineff.: ","Random noise:","Transient ineff.: "),
    Distribution = c("normal", "half-normal  ","normal", "half-normal  "),
    Assumption = Assumptions
  )
  print(a1[obj$toinclude,], quote = FALSE, right = F)

  # cat("\n------------")
  # cat(" Summary of the panel data: ------------\n\n", sep = "")
  # cat("\n------------ Summary of the panel data: ----------\n\n", sep = "")
  est.spd.left <- floor( (max.name.length+42-29) / 2 )
  est.spd.right <- max.name.length+42-29 - est.spd.left
  # cat("\n------------")
  # cat(" Summary of the panel data: ------------\n\n", sep = "")
  # cat("\n------------ Summary of the panel data: ----------\n\n", sep = "")
  cat("\n",rep("-", est.spd.left)," Summary of the panel data: ",rep("-", est.spd.right),"\n\n", sep ="")
  cat("   Number of obs       (NT) =",obj$nt,"", "\n")
  cat("   Number of groups     (N) =",obj$n,"", "\n")
  cat("   Obs per group: (T_i) min =",obj$dat.descr[3],"", "\n")
  cat("                        avg =",obj$dat.descr[4],"", "\n")
  cat("                        max =",obj$dat.descr[5],"", "\n")

  # cat("\n---------------")
  # cat(" Estimation results: ----------------\n\n", sep = "")
  # .printgtresfhet(obj$results, digits = digits, obj$Kb, obj$Kv0, obj$Ku0, obj$Kvi, obj$Kui, na.print = "NA")
  est.rez.left <- floor( (max.name.length+42-22) / 2 )
  est.rez.right <- max.name.length+42-22 - est.rez.left
  cat("\n",rep("-", est.rez.left)," Estimation results: ",rep("-", est.rez.right),"\n\n", sep ="")
  # cat("\n--------------- Estimation results: --------------\n\n", sep = "")
  .printgtresfhet(obj$results, digits = digits, obj$Kb, obj$Kv0, obj$Ku0, obj$Kvi, obj$Kui, na.print = "NA", max.name.length)
  # cat("____________________________________________________\n")

  if(obj$prod){
    myeff <- "technical"
    ndots <- "...."
    e_i0 <- obj$te_i0
    e_it <- obj$te_it
    e_ov <- obj$te_over
  }
  else {
    myeff <- "cost"
    ndots <- "........."
    e_i0 <- obj$ce_i0
    e_it <- obj$ce_it
    e_ov <- obj$ce_over
  }

  # sum.te_i0 <- paste0(" Summary of persistent ",myeff," efficiencies: ")
  # est.cle.left <- floor( (max.name.length+42-nchar(sum.te_i0)-1) / 2 )
  # est.cle.right <- max.name.length+42-nchar(sum.te_i0) -1- est.cle.left
  # if(est.cle.left <= 0) est.cle.left <- 1
  # if(est.cle.right <= 0) est.cle.right <- 1
  # cat("\n",rep("-", est.cle.left),sum.te_i0,rep("-", est.cle.right),"\n\n", sep ="")
  # .su(e_i0, print = TRUE)
  #
  # sum.te_it <- paste0(" Summary of transient ",myeff," efficiencies: ")
  # est.cle.left <- floor( (max.name.length+42-nchar(sum.te_it)-1) / 2 )
  # est.cle.right <- max.name.length+42-nchar(sum.te_it) -1- est.cle.left
  # if(est.cle.left <= 0) est.cle.left <- 1
  # if(est.cle.right <= 0) est.cle.right <- 1
  # cat("\n",rep("-", est.cle.left),sum.te_it,rep("-", est.cle.right),"\n\n", sep ="")
  # .su(e_it, print = TRUE)
  #
  # sum.te_over <- paste0(" Summary of overall ",myeff," efficiencies: ")
  # est.cle.left <- floor( (max.name.length+42-nchar(sum.te_over)-1) / 2 )
  # est.cle.right <- max.name.length+42-nchar(sum.te_over) -1- est.cle.left
  # if(est.cle.left <= 0) est.cle.left <- 1
  # if(est.cle.right <= 0) est.cle.right <- 1
  # cat("\n",rep("-", est.cle.left),sum.te_over,rep("-", est.cle.right),"\n\n", sep ="")
  # .su(e_ov, print = TRUE)

  sum.te <- paste0(" Summary of ",myeff," efficiencies: ")
  est.cle.left <- floor( (max.name.length+42-nchar(sum.te)-1) / 2 )
  est.cle.right <- max.name.length+42-nchar(sum.te) -1- est.cle.left
  if(est.cle.left <= 0) est.cle.left <- 1
  if(est.cle.right <= 0) est.cle.right <- 1
  cat("\n",rep("-", est.cle.left),sum.te,rep("-", est.cle.right),"\n\n", sep ="")
  .su(list(e_i0,e_it,e_ov), print = TRUE, width = 5, format = "fg", drop0trailing = FALSE, names = c("Persistent", "Transient", "Overall"))

  # ME
  if(obj$Ku0 > 1){
    cat("",rep("_", max.name.length+42-1),"", "\n", sep = "")

        sum.me_i0 <- paste0(" Summary of ME of Zi0 on persistent ",myeff," eff-s: ")
    est.cle.left <- floor( (max.name.length+42-nchar(sum.me_i0)-1) / 2 )
    est.cle.right <- max.name.length+42-nchar(sum.me_i0) -1- est.cle.left
    if(est.cle.left <= 0) est.cle.left <- 2
    if(est.cle.right <= 0) est.cle.right <- 2
    cat("\n",rep("-", est.cle.left),sum.me_i0,rep("-", est.cle.right),"\n\n", sep ="")
    .su(obj$me_i0, print = TRUE)

    sum.me_i0 <- paste0(" Summary of elasticity of E(u_i0) w.r.t Zi0: ")
    est.cle.left <- floor( (max.name.length+42-nchar(sum.me_i0)-1) / 2 )
    est.cle.right <- max.name.length+42-nchar(sum.me_i0) -1- est.cle.left
    if(est.cle.left <= 0) est.cle.left <- 2
    if(est.cle.right <= 0) est.cle.right <- 2
    cat("\n",rep("-", est.cle.left),sum.me_i0,rep("-", est.cle.right),"\n\n", sep ="")
    .su(obj$me_i0_elast, print = TRUE)
  }
  if(obj$Kui > 1){
    cat("",rep("_", max.name.length+42-1),"", "\n", sep = "")

    sum.me_i0 <- paste0(" Summary of ME of Zit on transient ",myeff," eff-s: ")
    est.cle.left <- floor( (max.name.length+42-nchar(sum.me_i0)-1) / 2 )
    est.cle.right <- max.name.length+42-nchar(sum.me_i0) -1- est.cle.left
    if(est.cle.left <= 0) est.cle.left <- 2
    if(est.cle.right <= 0) est.cle.right <- 2
    cat("\n",rep("-", est.cle.left),sum.me_i0,rep("-", est.cle.right),"\n\n", sep ="")
    .su(obj$me_it, print = TRUE)

    sum.me_i0 <- paste0(" Summary of elasticity of E(u_it) w.r.t Zit: ")
    est.cle.left <- floor( (max.name.length+42-nchar(sum.me_i0)-1) / 2 )
    est.cle.right <- max.name.length+42-nchar(sum.me_i0) -1- est.cle.left
    if(est.cle.left <= 0) est.cle.left <- 2
    if(est.cle.right <= 0) est.cle.right <- 2
    cat("\n",rep("-", est.cle.left),sum.me_i0,rep("-", est.cle.right),"\n\n", sep ="")
    .su(obj$me_it_elast, print = TRUE)
  }

  invisible( obj )
}