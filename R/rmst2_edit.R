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