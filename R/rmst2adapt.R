#' @name rmst2adapt
#' @aliases rmst2adapt
#' @title Restricted Mean Survival Times with Adaptively Selected Truncation Time Point
#' @description This function performs the procedure proposed by Horiguchi et al. (2018) <doi:10.1002/sim.7661>.
#' The method specifies a set of truncation time points, taus, for calculating restricted mean survival times (RMST),
#' performs testing for equality of the RMSTs between the two groups, and estimates the difference in RMST between the two groups at the specified tau's.
#' The multiplicity as a result of specifying several taus is taken into account in this procedure.

#' @usage  rmst2adapt(indata, tau_star, method="perturbation", nmethod=100000,
#'         seed=NULL, test="2_side", conf.int=0.95)
#' @param indata A data matrix (data frame). The 1st column is the time-to-event variable, the 2nd column is the event indicator (1=event, 0=censor), and the 3rd column is the treatment indicator (1=treatment, 0=control).
#' No missing values are allowed in this data matrix.
#' @param tau_star A vector indicating a set of tau's. All elements in \code{tau_star} need to be shorter than or equal to the minimum time point of the largest observed time in each of the two groups.
#' @param method A type of the resampling method. It supports \code{"perturbation"} (default).
#' @param nmethod A number of iterations for resampling. It is recommended to specify at least 100000 (default) or larger.
#' @param seed An integer value used for random number generation in the resampling procedures. Default is \code{NULL}.
#' @param test Specify \code{"1_side"} for the one-sided test where the alternative hypothesis is that the treatment group is superior to the control group with respect to survival time.
#' Specify \code{"2_side"} for the two-sided test where the alternative hypothesis is that the treatment group is not equal to the control group with respect to survival time.
#' Default is \code{"2_side"}.
#' @param conf.int Specify a confidence coefficient for calculating confidence bands for the differences in RMST. Default is \code{0.95}.
#' @return an object of class rmst2adapt.
#' @return \item{method}{The resampling method used in the analyses}
#' @return \item{nmethod}{The number of iterations for the resampling}
#' @return \item{test}{The type of test used in the analyses}
#' @return \item{candidate_taus}{The set of taus used in the analyses}
#' @return \item{observed_z}{The observed test statistic Z_star}
#' @return \item{p_value}{The p-value of testing for equality of the RMSTs between the two groups}
#' @return \item{conf_band}{The difference in RMST between the two groups at the specified taus}
#' @return \item{selected_tau}{The value of tau selected to summarize the treatment effect}
#' @references Horiguchi M, Cronin A, Takeuchi M, Uno H. A flexible and coherent test/estimation procedure based on restricted mean survival times for censored time-to-event data
#' in randomized clinical trials. Statistics in Medicine 2018. doi:10.1002/sim.7661.
#' @author Miki Horiguchi, Hajime Uno
#' @examples
#' #--- sample data ---#
#' data    = rmst2adapt.sample.data()
#' nmethod = 100 #This is only for example use.
#'               #Recommended to specify at least 100000 (default) or larger.
#'
#' a = rmst2adapt(indata=data, tau_star=seq(6,12,2), method="perturbation",
#'                nmethod=nmethod, test="2_side", seed=123)
#' print(a)

#' @export
############################################################################################
# main function 1: rmst2adapt()
############################################################################################
rmst2adapt <- function(indata, tau_star, method="perturbation", nmethod=100000, seed=NULL, test="2_side", conf.int=0.95){
  ##--observed Z*
  a = rmst2_edit(indata$time, indata$status, indata$arm, seq_tau=tau_star, type="observation", test=test)
  obs_rmstdiff    = a$unadjusted.result[,1]
  obs_rmstdiff_se = a$unadjusted.result[,2]
  obs_rmst1       = a$unadjusted.result[,6]
  obs_rmst0       = a$unadjusted.result[,7]
  
  check_box = matrix(NA, length(tau_star), 2)
  colnames(check_box) = c("tau", "Z_tau")
  check_box[,1] = tau_star
  check_box[,2] = a$unadjusted.result[,4] #rmst.diff.z, if one_side. abs(rmst.diff.z), if two_side.
  
  max_tau = check_box[,1][check_box[,2]==max(check_box[,2])]
  if(length(max_tau)>1){
    max_tau = max_tau[length(max_tau)]
  }
  
  obs_z_star = check_box[,2][check_box[,1]==max_tau]
  
  if(!is.null(seed)) set.seed(seed)
  
  if(method=="perturbation"){
    ##===================================================
    # resampling (Perturbation)
    ##===================================================
    pert_z_star = NULL
    for (g in 1:nmethod){
      b = rmst2_edit(time=indata$time, status=indata$status, arm=indata$arm, seq_tau=tau_star, type="perturbation", test=test, obs_rmstdiff=obs_rmstdiff)
      pert_check_box = matrix(NA, length(tau_star), 2)
      colnames(pert_check_box) = c("tau", "Z_tau")
      pert_check_box[,1] = tau_star
      pert_check_box[,2] = b$unadjusted.result[,4] #rmst.diff.z, if one_side. abs(rmst.diff.z) if two_side.
      
      pert_z_star[g] = max(pert_check_box[,2])[1]
    }
    
    #get p-value
    pval = sum(pert_z_star > obs_z_star)/nmethod
    
    #get CI
    c_alpha = quantile(pert_z_star, conf.int)
    
    #-- Get confidence band --
    rmstdiff_cb_low = obs_rmstdiff - obs_rmstdiff_se*c_alpha
    rmstdiff_cb_upp = obs_rmstdiff + obs_rmstdiff_se*c_alpha
    
    if(test=="1_side"){
      rmstdiff_cb_upp = tau_star
    }
    
  }
  
  #results at selected tau
  tmp8 = matrix(0, length(tau_star), 6)
  colnames(tmp8) = c("Tau", "RMST(arm1)", "RMST(arm0)", "RMST(arm1-arm0)", paste0("lower ", conf.int), paste0("upper ", conf.int))
  
  tmp8[,1] = tau_star
  tmp8[,2] = obs_rmst1
  tmp8[,3] = obs_rmst0
  tmp8[,4] = obs_rmstdiff
  tmp8[,5] = rmstdiff_cb_low
  tmp8[,6] = rmstdiff_cb_upp
  
  #output
  Z = list()
  Z$method          = method
  Z$nmethod         = nmethod
  Z$test            = test
  Z$candidate_taus  = tau_star
  Z$observed_z      = obs_z_star
  Z$p_value         = pval
  Z$conf_band       = tmp8
  Z$selected_tau    = max_tau
  
  class(Z) = "rmst2adapt"
  
  Z
}