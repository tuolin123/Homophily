### fit models for simulated data
check.AN = function(M = 10, # Simulation size 
                    n = 50,
                    k_0 = 10,
                    k_tot = 10,
                    Beta00 = 0.5,
                    Beta10 = 0.5,
                    size1 = 300,
                    size2 = 1100,
                    midsize = 200,
                    
                    Beta0_matrix, 
                    Beta1_matrix) {
  
  Beta00 = Beta00
  Beta10 = Beta10
  
  Beta0_chisq.MN = Beta1_chisq.MN  = c()
  Beta_est.BN = Var_est.BN = list()
  Beta_est.MN = Var_est.MN = list()
  
  for (m in 1:M) {
    
    ### Generating Network ###
    SN1 = network.sim(n = n,
                      k_0 = k_0,
                      k_tot = k_tot,
                      Beta00 = Beta00,
                      Beta10 = Beta10,
                      size1 = size1,
                      size2 = size2,
                      midsize = midsize) 
    
    x_jt_minus_x_new = SN1$x_jt_minus_x_new # dim: [n-(k_t-k_0)] * (k_t) or (k_t-1)
    Z_vec = SN1$Z_vec
    
    XX = plyr::ldply(x_jt_minus_x_new, rbind) ## Covariate
    XX = as.data.frame(XX)
    XX = as.matrix(XX)
    
    ############## BN ###############
    ZZ.BN = c()
    t_remove = c()
    for (t in 1:length(Z_vec)){
      if(Z_vec[[t]][1] != 1){
        update = Z_vec[[t]][-1]
        length(update) = max(sapply(Z_vec, length)-1)
        ZZ.BN = rbind(ZZ.BN, update)
      }else{
        t_remove = c(t_remove, t)
      }
    }
    ZZ.BN = as.matrix(ZZ.BN)
    if(is.null(t_remove)){
      XX = XX
    }else{
      XX = XX[-t_remove, ]
    }
    
    ### Inf loop for BN
    Rho_k.est.mean =  as.numeric()
    Beta.est.BN = matrix(0, nrow=ncol(XX), ncol=2)
    for (k in (1:ncol(XX))) {
      datfit = data.frame(cbind(ZZ.BN[, k], XX[, k]))
      datfit = datfit[complete.cases(datfit),]
      colnames(datfit) = c("y", "x")
      MLE2 = glm(y ~ x, data = datfit, family = "binomial")
      Beta.est.BN[k, ] = MLE2$coef
      
      l.beta.est.BN = Beta.est.BN[k,1] + Beta.est.BN[k,2] * XX[, k]
      Rho_k.est = exp(l.beta.est.BN)/(1+exp(l.beta.est.BN)) # BERNOULLI RHO
      Rho_k.est.mean[k] = mean(Rho_k.est, na.rm = T)
      
    }

    # BN for ALL:
    datfit = data.frame(cbind(c(ZZ.BN), c(XX)))
    datfit = datfit[complete.cases(datfit),]
    colnames(datfit) = c("y", "x")
    
    MLE2.all = glm(y ~ x, data = datfit, family = "binomial")
    Beta.est.BN.all = MLE2.all$coef
    Rho.est.BN = as.numeric()
    for(k in (1:ncol(XX))){
      l.beta.est.BN.all = Beta.est.BN.all[1] + Beta.est.BN.all[2] * XX[, k]
      Rho.est.BN[k] = mean(exp(l.beta.est.BN.all)/(1+exp(l.beta.est.BN.all)), na.rm = T) # BERNOULLI RHO
    }
    
    #BN.p_k.est = round(colMeans(Rho.est.BN, na.rm = T),4)
    BN.p_k.est = round(Rho.est.BN, 4)
    
    Beta_est.BN[[m]] = as.vector(Beta.est.BN.all)
    
    
    ############## MN ###############
    XX = plyr::ldply(x_jt_minus_x_new, rbind) ## Cov
    XX = as.data.frame(XX)
    XX = as.matrix(XX)
    
    ZZ.MN = c()
    for (t in 1:length(Z_vec)){
      update = as.vector(Z_vec[[t]])
      length(update) = max(sapply(Z_vec, length))
      ZZ.MN = rbind(ZZ.MN, update)
    }
    ZZ.MN = as.matrix(ZZ.MN)
    
    ### Inf loop of MN fit ###
    beta = c(0.5, 0.5)
    vec_len = sapply(Z_vec, length)
    MLE1.MN = optim(beta, logit.nll.all, XX=XX, ZZ=ZZ.MN,
                    len = vec_len, ind = 6:11,
                    method = "BFGS", hessian=TRUE)
    
    
    Beta.est.MN = MLE1.MN$par
    
    ### MN proportions ###
    l.beta.est.MN = Beta.est.MN[1] + Beta.est.MN[2] * XX
    exp.l.beta.MN = exp(l.beta.est.MN)
    p_k.est =  matrix(as.numeric(0), nrow=nrow(XX), ncol=ncol(XX))
    for(i in (1:nrow(XX))){
      col_k = length(exp.l.beta.MN[i,])
      p_k.est[i,] = exp.l.beta.MN[i,]/(1+sum(exp.l.beta.MN[i,],na.rm = T))
      
    }
    MN.p_k.est = round(colMeans(p_k.est, na.rm=T),4)
    
    Beta_est.MN[[m]] = MLE1.MN$par
  }
  
  Beta.emp.BN = rowMeans(as.data.frame(Beta_est.BN))
  Beta.emp.MN = rowMeans(as.data.frame(Beta_est.MN))
  
  
  return(list( 
    
    ### BN ###
    Rho_k.est.mean = Rho_k.est.mean, # BN for each column
    BN.p_k.est = BN.p_k.est,         # BN for ALL
    Beta.emp.BN = Beta.emp.BN,
    
    ### MN ###
    MN.p_k.est = MN.p_k.est,
    Beta.emp.MN = Beta.emp.MN
  ))
}
