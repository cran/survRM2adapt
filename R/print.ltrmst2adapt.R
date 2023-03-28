#' @name print.ltrmst2adapt
#' @aliases print.ltrmst2adapt
#' @title print.ltrmst2adapt
#' @description S3 method for class 'ltrmst2adapt'
#' @param x Object to be printed
#' @param digits Integer indicating the number of decimal places
#' @param ... Further arguments ignored in this function
#' @return returns summary output for class 'ltrmst2adapt'
#' @export
######################################
# print.ltrmst2adapt (hidden)
######################################
print.ltrmst2adapt <- function(x, digits=3, ...){
  
  cat("\n")
  
  diff10        = x$diff10_selected
  selected_tau1 = diff10["tau1"]
  tau2          = diff10["tau2"]
  
  #---
  arm0 = x$arm0[x$arm0[,"tau1"]==selected_tau1,"RMST"]
  arm1 = x$arm1[x$arm1[,"tau1"]==selected_tau1,"RMST"]
  names(arm0) = "RMST(arm0)"
  names(arm1) = "RMST(arm1)"
  
  cat("RMST over the time window [", selected_tau1,  ",", tau2, "] \n", sep="")
  cat("\n")
  rmst         = round(c(arm1, arm0), digits=digits)
  rmst.diff    = round(diff10[c(3,6,7,8)], digits=digits)
  
  
  print(rmst)
  cat("\n")
  print(rmst.diff)
  
  cat("\n\n")
  invisible(x)
}