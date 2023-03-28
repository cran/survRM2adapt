#=============================================
# twoprmst.R
# by Miki Horiguchi, Lu Tian,and Hajime Uno
#=============================================
twoprmst <- function(y, d, trt, tau1, tau2=NULL, B=50000, direction="from tau", seed=NULL, test="2_side", conf.int=0.95){
  #================
  #-- tau check --
  #================
  if(!is.null(tau2) & max(tau1)>=tau2){
        stop("tau1 needs to be shorter than tau2")
  }
 
  #=====================
  #-- max(tau2) check --
  #=====================
  idx = trt==0; yy = y[idx]; yy0max = max(yy); ss = d[idx]; ss0max = min(ss[yy==yy0max]);
  idx = trt==1; yy = y[idx]; yy1max = max(yy); ss = d[idx]; ss1max = min(ss[yy==yy1max]);

  yymax = max(yy0max, yy1max)
  yymin = min(yy0max, yy1max)
  
  #--Case 1: the last obs (smaller one)=event, the last obs (longer one)=event
  if(ss0max==1 & ss1max==1){
    if(!is.null(tau2)){
      if(tau2>yymax){stop(paste("The truncation time, tau2, needs to be shorter than or equal to ", round(yymax, digits=2)))}
      if(tau2<=yymax){NOTE=paste("The truncation time: tau2 =", tau2, " was specified.")}
    }else{
      tau2 = yymax
      NOTE = (paste("The truncation time, tau2, was not specified. Thus, the default tau2 ", round(yymax, digits=2)," is used."))
    }
  }
  
  #--Case 2: the last obs (smaller one)=event, the last obs (longer one)=censor
  if((ss0max==0 & ss1max==1 & yy0max>=yy1max) | (ss0max==1 & ss1max==0 & yy1max>yy0max)){
    if(!is.null(tau2)){
      if(tau2>yymax){stop(paste("The truncation time, tau2, needs to be shorter than or equal to ", round(yymax, digits=2)))}
      if(tau2<=yymax){NOTE=paste("The truncation time: tau2 =", tau2, " was specified.")}
    }else{
      tau2 = yymax
      NOTE = paste("The truncation time, tau2, was not specified. Thus, the default tau2 ", round(yymax, digits=2)," is used.")
    }
  }
  
  #--Case 3: the last obs (smaller one)=censor, the last obs (longer one)=event
  if((ss0max==1 & ss1max==0 & yy0max>=yy1max) | (ss0max==0 & ss1max==1 & yy1max>yy0max)){
    if(!is.null(tau2)){
      if(tau2>yymin){stop(paste("The truncation time, tau2, needs to be shorter than or equal to ", round(yymin, digits=2)))}
      if(tau2<=yymin){NOTE=paste("The truncation time: tau2 =", tau2, " was specified.")}
    }else{
      tau2 = yymin
      NOTE = (paste("The truncation time, tau2, was not specified. Thus, the default tau2 ", round(yymin, digits=2)," is used."))
    }
  }
  
  #--Case 4: the last obs (smaller one)=censor, the last obs (longer one)=censor
  if(ss0max==0 & ss1max==0){
    if(!is.null(tau2)){
      if(tau2<=yymin){
        NOTE=paste("The truncation time: tau2 =", tau2, " was specified.")
      }
      if(tau2>yymin){
        stop(paste("The truncation time, tau2, needs to be shorter than or equal to the minimum of the largest observed time on each of the two groups: ", round(yymin, digits=2)))
      }
    }else{
      tau2 = yymin
      NOTE = (paste("The truncation time, tau2, was not specified. Thus, the default tau2 ", round(yymin, digits=2)," is used."))
    }
  }
  
  #==========
  #-- main --
  #==========
  fit1  = prmst_mod(y[trt==1], d[trt==1], tau1=tau1, tau2=tau2, direction=direction)
  fit0  = prmst_mod(y[trt==0], d[trt==0], tau1=tau1, tau2=tau2, direction=direction)
  delta = fit1$auc-fit0$auc
  sigma = fit1$sigma+fit0$sigma
  se    = sqrt(diag(sigma))

  if(!is.null(seed)){
  	set.seed(seed)
  }
  n         = length(delta)
  cor.sigma = sigma
  for(i in 1:n){
    for(j in 1:n){
      cor.sigma[i,j]=cor.sigma[i,j]/prod(se[c(i,j)])
    }
  }

  error = rmvnorm(n=B, mean=rep(0,n), sigma=cor.sigma, method="eigen")
  
  #----------
  alpha = 1-conf.int
  if(test=="2_side"){
  	cut = quantile(apply(abs(error), 1, max), 1-alpha)
  	pv  = mean(apply(abs(error), 1, max) > max(abs(delta)/se))
    cr  = cbind(delta-cut*se, delta+cut*se)
  }else{
  	#--If test=="1_side"
  	cut = quantile(apply(error,1,max), 1-alpha)
  	pv  = mean(apply(error, 1, max) > max(delta/se))
    cr  = cbind(delta-cut*se, fit1$tau2)
  }
  
  #=============
  #-- output1 --
  #=============
  Z = list()
  Z$tau1        = fit1$tau1
  Z$tau2        = fit1$tau2
  Z$direction   = direction
  Z$note        = NOTE
  Z$result_arm1 = fit1$out
  Z$result_arm0 = fit0$out
  Z$cut         = cut
  Z$delta       = delta
  Z$se          = se
  Z$sigma       = sigma
  Z$cr          = cr
  Z$pv          = pv
  if(test=="2_side"){
  	out = cbind(Z$tau1, Z$tau2, delta, se, abs(delta/se), cr)
  }else{
  	out = cbind(Z$tau1, Z$tau2, delta, se, delta/se, cr)
  }
  colnames(out) = c("tau1", "tau2", "RMST(arm1-arm0)", "SE", "Z", paste0(c("lower ", "upper "), (1-alpha)))
  
  Z$out         = out
  Z$chosen      = out[out[,"Z"]==max(out[,"Z"]),]
  Z$seed        = seed
  Z$test        = test
  Z$conf.int    = conf.int
  return(Z)
}

