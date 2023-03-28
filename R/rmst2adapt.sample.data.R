#' @name rmst2adapt.sample.data
#' @aliases  rmst2adapt.sample.data
#' @title Generate a sample data from the pbc data
#' @description This is a function to retrieve 312 randomized patients from the pbc data in survival package.
#' @usage rmst2adapt.sample.data(t.unit="year")
#' @param t.unit Specify the time unit. It supports "year" (default), "month", and "day".
#' @return returns a data frame
#' @details The function creates a sample dataset to illustrate the usage of the function \code{rmst2adapt()} in this package.
#' The original pbc data in \code{survival} package consists of 418 patients data.
#' This function loads the pbc data, select the 312 patients who were randomized.
#' The status variable is edited, so that 1 indicates death and 0 indicates alive.
#' @seealso \code{pbc} in survival package
#' @examples
#' D=rmst2adapt.sample.data()
#' head(D)
#' @export
#######################################
# rmst2adapt sample data
#######################################
rmst2adapt.sample.data <- function(t.unit="year"){
  tmp = survival::pbc
  D   = tmp[1:312,c(2:4)]
  
  if(t.unit=="year"){
    D$time=D$time/365.25
  }
  if(t.unit=="month"){
    D$time=D$time/365.25*12
  }
  if(t.unit=="day"){
    D$time=D$time
  }
  
  D$status = as.numeric(D$status==2)
  D$arm    = as.numeric(D$trt==1)
  
  DA = D[,-3]
  DA
}
NULL