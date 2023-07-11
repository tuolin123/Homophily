# function for multinomial regression MLE  

logit.nll = function (beta, XX, ZZ) {
  l.beta.est = beta[1] + beta[2]*XX 
  upper.est = exp(l.beta.est)
  Denomenator = as.matrix(1+rowSums(upper.est, na.rm = T)) #Sum_(h=1)^3 exp(X * Beta_(h))    

  #LL=0
  Pr_1 = 1/Denomenator
  L_1 = colSums(ZZ[,1] * log(Pr_1), na.rm = T)
  LL=0+L_1
  for (j in 1:ncol(upper.est)){
    Pr_j = upper.est[,j]/Denomenator
    L_j = colSums(ZZ[,j+1] * log(Pr_j), na.rm = T)
    LL = LL + L_j #+ L_1
  }
  return(-LL)
  
}

logit.nll.all = function(beta, XX, ZZ, len, ind){
  ll = 0
  for(k in ind){
    X_comp = XX[len == k,]
    X_comp = X_comp[,colSums(is.na(X_comp))<nrow(X_comp)]
    ZZ_comp = ZZ[len == k,]
    ZZ_comp = ZZ_comp[,colSums(is.na(ZZ_comp))<nrow(ZZ_comp)]
    ll1 = logit.nll(beta, X_comp, ZZ_comp)
    ll = ll + ll1
  }
  return(ll)
}
