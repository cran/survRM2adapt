#=============================================
# prmst_mod.R
# by Miki Horiguchi, Lu Tian,and Hajime Uno
# last updated: 2019/06/21
#=============================================

#-------------------
prmst_mod <- function(y, d, tau1, tau2, direction="from tau"){
  
  fit = survfit(Surv(y, d)~1)
  tt  = c(0, summary(fit)$time)
  sp  = c(1, summary(fit)$surv)
  
  m   = length(tt)
  n   = length(tau1)
  
  #--- RMST from 0 to tau2 (est)--
  index         = (1:m)[tt<=tau2]
  tt.interval   = c(0, tt[index], tau2)
  prob.interval = c(1, sp[index])
  auc_tau2      = sum(diff(tt.interval)*prob.interval)
  
  #--- RMST from 0 to each of the tau1's and tau2 (est) 
  #--- RMST from each of the tau1's to tau2 (est)--
  auc_tau1 = rep(0, n+1)
  auc      = rep(0, n)
  for(i in 1:n){
    index         = (1:m)[tt<=tau1[i]]
    tt.interval   = c(0, tt[index], tau1[i])
    prob.interval = c(1, sp[index])
    auc_tau1[i]   = sum(diff(tt.interval)*prob.interval)
    auc[i]        = auc_tau2 - auc_tau1[i]
  }
  auc_tau1[n+1]=auc_tau2
  
  #--- RMST from 0 to tau2 (var)--
  int2 = rep(0, m)
  for(b in 1:m){
    index         = (1:m)[tt>=tt[b] & tt<=tau2]
    tt.interval   = c(tt[index], tau2)
    prob.interval = sp[index]
    int2[b]       = sum(diff(tt.interval)*prob.interval)
  }
  
  nevent = nrisk = rep(0, m)
  for(b in 1:m){
    nevent[b] = sum(y[d==1]==tt[b])
    nrisk[b]  = sum(y>=tt[b])
  }
  
  idv      = int2^2*nevent/nrisk^2
  var_tau2 = sum(idv[tt<tau2])
  
  #--- RMST from 0 to tau1 and tau2 (variance-covariance matrix)--
  tau12    = c(tau1,tau2)
  int1_mat = matrix(0, nrow=n+1, ncol=m)
  for(i in 1:(n+1)){
    for(b in 1:m){
      index         = (1:m)[tt>=tt[b] & tt<=tau12[i]]
      tt.interval   = c(tt[index], tau12[i])
      prob.interval = sp[index]
      int1_mat[i,b] = sum(diff(tt.interval)*prob.interval)
    }
  }
  
  sigma_tau12   = matrix(0, nrow=n+1, ncol=n+1)
  cov_tau1_tau2 = rep(0,n)
  for(i in 1:(n+1)){
    for(j in 1:(n+1)){
      idv              = int1_mat[i,]*int1_mat[j,]*nevent/nrisk^2
      sigma_tau12[i,j] = sum(idv[tt<tau12[min(i,j)]])
    }
  }
  
  #--- Adaptive RMST (Horiguchi et al. Statistics in Medicine (2018)) ---
 # if(direction=="from 0"){
 #   start_tau     = rep(0, n+1)
 #   end_tau       = c(tau1, tau2)
 #   auc           = auc_tau1
 #   se            = sqrt(diag(sigma_tau12))
 #   sigma         = sigma_tau12
 #   out           = cbind(start_tau, end_tau, auc, se)
 #   colnames(out) = c("tau1", "tau2", "RMST", "SE")
 #} 
  
  #--- Adaptive long-term RMST (Horiguchi, Tian, and Uno. Statistics in Medicine (2023) ---
  if(direction=="from tau"){
    start_tau = tau1
    end_tau   = rep(tau2, n)
    auc       = auc
    
    #--- RMST from tau1 to tau2 (var-cov matrix)--
    sigma = matrix(0,n,n)
    for(i in 1:n){
      for(j in 1:n){
        sigma[i,j] = var_tau2 + sigma_tau12[i,j] -  sigma_tau12[n+1,i] - sigma_tau12[n+1,j]
      }
    }
    se            = sqrt(diag(sigma))
    out           = cbind(start_tau, end_tau, auc, se)
    colnames(out) = c("tau1","tau2","RMST","SE")
  }
  
  #==================    
  Z = list()
  Z$tau1      = start_tau
  Z$tau2      = end_tau
  Z$auc       = auc
  Z$se        = se
  Z$sigma     = sigma
  Z$out       = out
  Z$direction = direction
  return(Z)
}
