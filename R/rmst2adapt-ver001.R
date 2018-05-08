#' @name survRM2adapt-package
#' @aliases  survRM2adapt-package
#' @docType  package
#' @title Flexible and Coherent Test/Estimation Procedure Based on Restricted Mean Survival Times
#' @description
#' Performs the procedure proposed by Horiguchi et al. (2018) <doi:10.1002/sim.7661>.
#' The method specifies a set of truncation time points tau's for calculating restricted mean survival times (RMST),
#' performs testing for equality, and estimates the difference in RMST between two groups at the specified tau's.
#' Multiplicity by specifying several tau's is taken into account in this procedure.
#'
#' @author Miki Horiguchi, Hajime Uno
#'
#' Maintainer: Miki Horiguchi <horiguchimiki@gmail.com>
#'
#' @references
#' Horiguchi M, Cronin A, Takeuchi M, Uno H. A flexible and coherent test/estimation procedure based on restricted mean survival times for censored time-to-event data
#' in randomized clinical trials. Statistics in Medicine 2018. doi:10.1002/sim.7661.
#'
#' @keywords
#' survival
#' @seealso
#' survival survRM2
#' @import survival
#' @importFrom stats pnorm qnorm rexp quantile
#' @importFrom utils data
#' @examples
#' #--- sample data ---#
#' data    = rmst2adapt.sample.data()
#' nmethod = 100 #This is only for example use.
#'               #Recommended to specify at least 100000 (default) or larger.
#'
#' a = rmst2adapt(indata=data, tau_star=seq(6,12,2), method="perturbation",
#'                nmethod=nmethod, test="2_side")
#' print(a)
NULL


#' @name rmst2adapt.sample.data
#' @aliases  rmst2adapt.sample.data
#' @title Generate a sample data from the pbc data
#' @description This is a function to retrieve 312 randomized patients from the pbc data in survival package.
#' @usage rmst2adapt.sample.data()
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
rmst2adapt.sample.data <- function(){
  tmp = survival::pbc
  D   = tmp[1:312,c(2:4)]

  D$time   = D$time/365.25
  D$status = as.numeric(D$status==2)
  D$arm    = as.numeric(D$trt==1)

  DA = D[,-3]
  DA
}
NULL





#' @name rmst2adapt
#' @aliases rmst2adapt
#' @title Flexible and Coherent Test/Estimation Procedure Based on Restricted Mean Survival Times
#' @description Performs the procedure proposed by Horiguchi et al. (2018) <doi:10.1002/sim.7661>.
#' The method specifies a set of truncation time points tau's for calculating restricted mean survival times (RMST),
#' performs testing for equality, and estimates the difference in RMST between two groups at the specified tau's.
#' Multiplicity by specifying several tau's is taken into account in this procedure.

#' @usage  rmst2adapt(indata, tau_star, method="perturbation", nmethod=100000,
#'         seed=NULL, test="2_side", conf.int=0.95)
#' @param indata A data matrix (data frame). The 1st column is time-to-event variable, the 2nd column is event indicator (1=event, 0=censor), and the 3rd column is the treatment indicator (1=treatment, 0=control).
#' No missing values are allowed in this data matrix.
#' @param tau_star A vector indicating a set of tau's. All elements in \code{tau_star} need to be shorter than or equal to the minimum of the largest observed time on each of the two groups.
#' @param method A name of the resampling method. It supports \code{"perturbation"} (default) and \code{"permutation"}.
#' @param nmethod A number of iterations for the resampling. Recommended to specify at least 100000 (default) or larger.
#' @param seed An integer value, used for the random number generation in the resampling procedures. Default is \code{NULL}.
#' @param test Specify \code{"1_side"} for the one-sided test where the alternative hypothesis is that treatment group is superior to control group with respect to survival.
#' Specify \code{"2_side"} for the two-sided test where the alternative hypothesis is that treatment group is not equal to control group with respect to survival.
#' Default is \code{"2_side"}.
#' @param conf.int Specify confidence coefficient for confidence bands of the differences in RMST. Default is \code{0.95}.
#' @return an object of class rmst2adapt.
#' @return \item{method}{The resampling method used in the analyses}
#' @return \item{nmethod}{The number of iterations for the resampling}
#' @return \item{test}{The type of test used in the analyses}
#' @return \item{candidate_taus}{The set of tau's used in the analyses}
#' @return \item{observed_z}{The observed test statistic Z_star}
#' @return \item{p_value}{The p-value of testing for equality}
#' @return \item{conf_band}{The difference in RMST between two groups at the specified tau's}
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
#'                nmethod=nmethod, test="2_side")
#' print(a)
NULL

############################################################################################
# rmst1_edit (one-arm) -- hidden
# Edited 'rmst1 by survRM2' to output the z value for both one_sided and two_sided test
############################################################################################
rmst1_edit <- function(time, status, seq_tau, weight=NULL){
  #-- time
  #-- statuts
  #-- seq_tau -- truncation times (vector)
  #-- weight=NULL --for perturbation

  ft = survfit(Surv(time, status)~1, weight=weight)

  rmst     = NULL
  rmst.var = NULL
  for(i in 1:length(seq_tau)){
    idx        = ft$time<=seq_tau[i]

    wk.time    = sort(c(ft$time[idx],seq_tau[i]))
    wk.surv    = ft$surv[idx]
    wk.n.risk  = ft$n.risk[idx]
    wk.n.event = ft$n.event[idx]

    time.diff  = diff(c(0, wk.time))
    areas      = time.diff * c(1, wk.surv)
    rmst[i]    = sum(areas)

    wk.var     =  ifelse((wk.n.risk-wk.n.event)==0, 0,
                         wk.n.event /(wk.n.risk *(wk.n.risk - wk.n.event)))
    wk.var     = c(wk.var,0)
    rmst.var[i] = sum(cumsum(rev(areas[-1]))^2 * rev(wk.var)[-1])
  }

  box = matrix(0,length(seq_tau),3)
  box[,1] = seq_tau
  box[,2] = rmst
  box[,3] = rmst.var
  colnames(box) = c("seq_tau", "rmst", "rmst.var")

  Z = list()
  Z$seq_tau  = box[,1]
  Z$rmst     = box[,2]
  Z$rmst.var = box[,3]
  return(Z)
}
NULL
############################################################################################
# rmst2_edit -- hidden
# Edited 'rmst2 by survRM2'
############################################################################################
rmst2_edit <- function(time, status, arm, seq_tau, type=NULL, test="2_side", obs_rmstdiff=NULL){
  #-- time
  #-- statuts
  #-- arm (1 or 0)
  #-- type ("perturbation", "permutation", or "observation")
  #-- test (1_side" or "2_side")
  #-- obs_rmstdiff --result of the difference in RMST via type="observation"

  #==================================
  #  initial check
  #==================================
  #--- tau ---
  idx=arm==0; tt=time[idx]; tau0max=max(tt)
  idx=arm==1; tt=time[idx]; tau1max=max(tt)
  tau_max = min(tau0max, tau1max)


  #---------------------
  if(!is.null(seq_tau)){
    if(seq_tau[length(seq_tau)] <= tau_max){
      NOTE=paste("A set of the truncation times: tau_star =", seq_tau, " was specified.")
    }
    if(seq_tau[length(seq_tau)]> tau_max){
      stop(paste("All elements in tau_star need to be shorter than or equal to the minimum of the largest observed time on each of the two groups: ", round(tau_max, digits=3)))
    }
  }


  if(type=="observation" | type=="permutation"){
    wk1 = rmst1_edit(time=time[arm==1], status=status[arm==1], seq_tau=seq_tau, weight=NULL)
    wk0 = rmst1_edit(time=time[arm==0], status=status[arm==0], seq_tau=seq_tau, weight=NULL)

    #--- contrast (RMST difference) ---
    rmst.diff.10     = wk1$rmst-wk0$rmst
    rmst.diff.10.se  = sqrt(wk1$rmst.var + wk0$rmst.var)
    rmst.diff.z      = rmst.diff.10/rmst.diff.10.se

    if(test=="1_side"){
      rmst.diff.z.1side      = rmst.diff.z
      rmst.diff.pval.1side   = 1-pnorm(rmst.diff.z.1side) # one-sided test (upper)
      rmst.diff.result.1side = cbind(rmst.diff.10, rmst.diff.10.se, 1, rmst.diff.z.1side, rmst.diff.pval.1side, wk1$rmst, wk0$rmst)
      #--- results ---
      out = rmst.diff.result.1side
    }else{
      #test=="2_side"
      rmst.diff.z.2side      = abs(rmst.diff.z)
      rmst.diff.pval.2side   = pnorm(-rmst.diff.z.2side)*2 # two-sided test
      rmst.diff.result.2side = cbind(rmst.diff.10, rmst.diff.10.se, 2, rmst.diff.z.2side, rmst.diff.pval.2side, wk1$rmst, wk0$rmst)
      #--- results ---
      out = rmst.diff.result.2side
    }
  }else{
    #type=="perturbation"
    n1  = length(time[arm==1])
    n0  = length(time[arm==0])
    wt1 = rexp(n1)
    wt0 = rexp(n0)
    wk1=rmst1_edit(time=time[arm==1], status=status[arm==1], seq_tau=seq_tau, weight=wt1)
    wk0=rmst1_edit(time=time[arm==0], status=status[arm==0], seq_tau=seq_tau, weight=wt0)

    #--- contrast (RMST difference) ---
    pert.rmst.diff.10    = (wk1$rmst-wk0$rmst) - obs_rmstdiff  #centering
    pert.rmst.diff.10.se = sqrt(wk1$rmst.var + wk0$rmst.var)
    pert.rmst.diff.z     = pert.rmst.diff.10/pert.rmst.diff.10.se

    if(test=="1_side"){
      pert.rmst.diff.z.1side      = pert.rmst.diff.z
      pert.rmst.diff.pval.1side   = 1-pnorm(pert.rmst.diff.z.1side) # one-sided test (upper)
      pert.rmst.diff.result.1side = cbind(pert.rmst.diff.10, pert.rmst.diff.10.se, 1, pert.rmst.diff.z.1side, pert.rmst.diff.pval.1side, wk1$rmst, wk0$rmst)
      #--- results ---
      out = pert.rmst.diff.result.1side
    }else{
      #test=="2side"
      pert.rmst.diff.z.2side      = abs(pert.rmst.diff.z)
      pert.rmst.diff.pval.2side   = pnorm(-pert.rmst.diff.z.2side)*2 # two-sided test
      pert.rmst.diff.result.2side = cbind(pert.rmst.diff.10, pert.rmst.diff.10.se, 2, pert.rmst.diff.z.2side, pert.rmst.diff.pval.2side, wk1$rmst, wk0$rmst)
      #--- results ---
      out = pert.rmst.diff.result.2side
    }
  }

  #--- results ---
  rownames(out) = paste0("tau",seq_tau)
  colnames(out) = c("Est.", "S.E.", "test-side", "z", "p", "rmst1", "rmst0")

  #--- output ---
  Z=list()
  Z$unadjusted.result = out

  Z
}
NULL

######################################
# shuffle.flat (ver.2) -- hidden
######################################
shuffle.flat <- function(data, key.var, seed=NULL, by_miss_pattern=FALSE){

  if(!is.null(seed)) set.seed(seed)

  #--- figure out the patterns ---
  tmp = data[,-1]
  n   = nrow(tmp)
  k   = ncol(tmp)
  pattern = rep(0, n)
  for (i in 1:k){
    pattern = pattern + as.numeric(!is.na(tmp[,i]))*(10^(i-1))
  }
  pattern

  unique_pattern = sort(unique(pattern))
  npatterns      = length(unique_pattern)

  tmp2 = cbind(data, pattern)
  tmp2

  #--- permuate by patterns ---
  if(by_miss_pattern==TRUE){
    D=c()
    for (i in 1:npatterns){
      idx = pattern == unique_pattern[i]
      tmp_D   = data[idx,]
      tmp_n   = nrow(tmp_D)
      key.org = tmp_D[,key.var]
      key.new = sample(key.org, size=tmp_n)
      tmp_D[,key.var] = key.new
      D = rbind(D, tmp_D)
    }
  }
  #--- permuate entire (ver.1)---
  if(by_miss_pattern==FALSE){
    D = data
    n = nrow(D)
    key.org = D[,key.var]
    key.new = sample(key.org, size=n)
    D[,key.var] = key.new
  }
  D
}
NULL

#' @export
############################################################################################
# main function: rmst2adapt
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

  }else{
    ##===================================================
    # resampling (Permutation)
    ##===================================================
    perm_z_star = NULL
    for (k in 1:nmethod){
      perm_dat = shuffle.flat(indata, "arm")

      b = rmst2_edit(time=perm_dat$time, status=perm_dat$status, arm=perm_dat$arm, seq_tau=tau_star, type="permutation", test=test)
      perm_check_box = matrix(NA, length(tau_star), 2)
      colnames(perm_check_box) = c("tau", "Z_tau")
      perm_check_box[,1] = tau_star
      perm_check_box[,2] = b$unadjusted.result[,4] #rmst.diff.z, if one_side. abs(rmst.diff.z) if two_side.

      c = max(perm_check_box[,2])
      if(length(c)>1){
        perm_z_star[k] = c[1]
      }else{
        perm_z_star[k] = c
      }
    }

    #get p-value
    pval = sum(perm_z_star > obs_z_star)/nmethod

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
NULL


#' @name print.rmst2adapt
#' @aliases print.rmst2adapt
#' @title print.rmst2adapt
#' @description S3 method for class 'rmst2adapt'
#' @param x Object to be printed
#' @param digits Integer indicating the number of decimal places
#' @param ... Further arguments ignored in this function
#' @export
######################################
# print.rmst2adapt (hidden)
######################################
print.rmst2adapt <- function(x, digits=3, ...){

  cat("\n")

  taus = x$candidate_taus
  pval = round(x$p_value, digits=digits)

  cat ("<Test result> \n")

  cat("Candidate tau's:", taus, "\n")
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
