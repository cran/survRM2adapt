#' @name print.rmst2adapt
#' @aliases print.rmst2adapt
#' @title print.rmst2adapt
#' @description S3 method for class 'rmst2adapt'
#' @param x Object to be printed
#' @param digits Integer indicating the number of decimal places
#' @param ... Further arguments ignored in this function
#' @return returns summary output for class 'rmst2adapt'
#' @export
######################################
# print.rmst2adapt (hidden)
######################################
print.rmst2adapt <- function(x, digits=3, ...){
  
  cat("\n")
  
  taus = x$candidate_taus
  pval = round(x$p_value, digits=digits)
  
  cat ("<Test result> \n")
  
  cat("Candidate values of tau1:", taus, "\n")
  cat("\n")
  cat("P-value:", pval,  "\n")
  
  cat("\n")
  cat("\n")
  
  cat ("<Treatment effect estimation> \n")
  
  rmst         = round(x$conf_band[x$conf_band[,1]==x$selected_tau,][c(2:3)], digits=digits)
  rmst.diff    = round(x$conf_band[x$conf_band[,1]==x$selected_tau,][c(4:6)], digits=digits)
  
  cat("Selected tau:", x$selected_tau, "\n")
  cat("\n")
  print(rmst)
  cat("\n")
  print(rmst.diff)
  
  cat("\n\n")
  invisible(x)
}
NULL