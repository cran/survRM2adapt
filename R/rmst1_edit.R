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