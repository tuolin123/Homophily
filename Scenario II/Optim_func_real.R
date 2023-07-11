# function for multinomial regression MLE  

logit.nll = function (beta, XX1, XX2, ZZ) {
  
  l.beta.est = beta[1] + beta[2]*XX1 + beta[3]*XX2
  upper.est = exp(l.beta.est)
  Denomenator = as.matrix(1+rowSums(upper.est, na.rm = T))  
  
  #LL=0
  Pr_1 = 1/Denomenator
  L_1 = colSums(ZZ[,1] * log(Pr_1), na.rm = T)
  LL=0+L_1
  for (j in 1:ncol(upper.est)){
    Pr_j = upper.est[,j]/Denomenator
    L_j = colSums(ZZ[,j+1] * log(Pr_j), na.rm = T)
    LL = LL + L_j
  }
  return(-LL)
}

logit.nll.all = function(beta, XX1, XX2, ZZ, len, ind){
  ll = 0
  for(k in ind){
    X1_comp = XX1[len == k,]
    X2_comp = XX2[len == k,]
    X1_comp = X1_comp[,colSums(is.na(X1_comp))<nrow(X1_comp)]
    X2_comp = X2_comp[,colSums(is.na(X2_comp))<nrow(X2_comp)]
    ZZ_comp = ZZ[len == k,]
    ZZ_comp = ZZ_comp[,colSums(is.na(ZZ_comp))<nrow(ZZ_comp)]
    ll1 = logit.nll(beta, X1_comp, X2_comp, ZZ_comp)
    ll = ll + ll1
  }
  return(ll)
}
