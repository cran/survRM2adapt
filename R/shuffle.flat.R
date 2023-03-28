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