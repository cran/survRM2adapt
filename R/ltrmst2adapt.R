#' @name ltrmst2adapt
#' @aliases ltrmst2adapt
#' @title Long-Term Restricted Mean Survival Times with Adaptively Selected Lower End of Time Window 
#' @description This function performs the method proposed by Horiguchi, Tian, and Uno. (2023) <doi:10.1002/sim.9662>.
#' Specifically, it estimates the restricted mean survival time (RMST) within a prespecified time window (from tau1 to tau2) to quantify the long-term treatment benefit. 
#' Instead of choosing one specific time point as the lower end of the time window (tau1), the procedure allows users to prespecify a set of time points.
#' The procedure picks one time point among the set of tau1 values that gives the most significant difference in the long-term RMST between the two groups.
#' It then performs testing for equality of the long-term RMSTs between the two groups and estimates the difference in RMST within the adaptively selected time window [tau1, tau2] between the two groups.
#' Multiplicity as a result of specifying several values for tau1 is taken into account in this procedure.


#' @usage  ltrmst2adapt(indata, tau1, tau2, iteration=50000, seed=NULL, test="2_side", conf.int=0.95)
#' @param indata A data matrix (data frame). The 1st column is the time-to-event variable, the 2nd column is the event indicator (1=event, 0=censor), and the 3rd column is the treatment indicator (1=treatment, 0=control).
#' No missing values are allowed in this data matrix.
#' @param tau1 An integer value or a vector indicating a set of tau1 values for the lower end of the time window.
#' @param tau2 An integer value indicating the upper end of the time window. When tau2 = NULL, the default value is used. See Details for the definition of the default value for tau2.
#' @param iteration A number of iterations for the resampling (the multivariate normal distribution-based perturbation method). It is recommended to specify at least 50000 (default) or larger.
#' @param seed An integer value used for random number generation in the resampling procedure. Default is \code{NULL}.
#' @param test Specify \code{"1_side"} for the one-sided test where the alternative hypothesis is that the treatment group is superior to the control group with respect to survival time. 
#' Specify \code{"2_side"} for the two-sided test where the alternative hypothesis is that the treatment group is not equal to the control group with respect to survival time.
#' Default is \code{"2_side"}.
#' @param conf.int Specify a confidence coefficient for calculating confidence bands for the differences in long-term RMST. Default is \code{0.95}.
#' @return an object of class ltrmst2adapt.
#' @return \item{iteration}{The number of iterations for resampling from a multivariate normal distribution with a mean 0 and a variance-covariance matrix.}
#' @return \item{test}{The type of test used in the analyses}
#' @return \item{arm1}{The RMST [tau1, tau2] estimation for arm1}
#' @return \item{arm0}{The RMST [tau1, tau2] estimation for arm0}
#' @return \item{diff10}{The difference in RMST [tau1, tau2] between the two groups (arm1 minus arm0)}
#' @return \item{diff10_selected}{The difference in RMST within the selected time window [tau1, tau2] between the two groups (arm1 minus arm0) with the normal confidence interval considering the randomness of selecting one time window. The p-value for the RMST difference test within the selected time window [tau1, tau2] is also provided.}
#' @references Horiguchi M, Tian L, Uno H. On assessing survival benefit of immunotherapy using long-term restricted mean survival time. Statistics in Medicine 2023. DOI:10.1002/sim.9662.
#' @author Miki Horiguchi, Lu Tian, Hajime Uno
#' @details The definition of the default value for tau2. Let x1 and x0 be the maximum observed time in Group 1 and Group 0, respectively. 
#' Case 1: If the last observations in Group 1 and Group 0 are "event," then tau = max(x1, x0). 
#' Case 2-1: If the last observation in Group 1 is "event," the last observation in Group 0 is "censor," and x1 <= x0, tau2 = max(x1, x0) = x0. 
#' Case 2-2: If the last observation in Group 0 is "event," the last observation in Group 1 is "censor," and x1 > x0, tau2 = max(x1, x0) = x1. 
#' Case 3-1: If the last observation in Group 1 is "event," the last observation in Group 0 is "censor," and x1 > x0, tau2 = min(x1, x0) = x0. 
#' Case 3-2: If the last observation in Group 0 is "event," the last observation in Group 1 is "censor," and x1 <= x0, tau2 = min(x1, x0) = x1. 
#' Case 4: If the last observations in Group 1 and Group 0 are "censor," then tau = min(x1, x0).
#' @examples
#' #--- sample data ---#
#' data = cm214_pfs
#' b    = ltrmst2adapt(indata=data, tau1=c(0,1,2,3), tau2=10, 
#'                test="2_side", seed=123)
#' print(b)


#' @export
ltrmst2adapt <- function(indata, tau1, tau2=NULL, iteration=50000, seed=NULL, test="2_side", conf.int=0.95){
  
  #-- Check tau1 and tau2
  if(length(tau2)>1){
    stop("This option is not available with this version.")
  }
  
  if(length(tau1)>1){
    method_type    = "RMST[tau1,tau2] with adaptively selected tau1"
    direction_type = "from tau"
  }
  
  if(length(tau1)==1){
    method_type    = "RMST[tau1,tau2] with fixed tau1"
    direction_type = "from tau"
  }

  
 
  #--adaptive LT-RMST
  wk = twoprmst(y=indata$time, d=indata$status, trt=indata$arm, tau1=tau1, tau2=tau2, B=iteration, direction=direction_type, seed=seed, test=test, conf.int=conf.int)
  
 
  #============
  #-- output --
  #============
  tmpout = c(wk$chosen, p_value=wk$pv)

  #---
  Z2 = list()
  #Z2$type   = method_type
  Z2$iteration       = iteration
  Z2$test            = test
  Z2$candidate_tau1  = tau1
  Z2$arm1            = wk$result_arm1
  Z2$arm0            = wk$result_arm0
  Z2$diff10          = wk$out[,c("tau1", "tau2", "RMST(arm1-arm0)", "SE", "Z")]
  Z2$diff10_selected = tmpout
  
  class(Z2) = "ltrmst2adapt"
  
  Z2
}