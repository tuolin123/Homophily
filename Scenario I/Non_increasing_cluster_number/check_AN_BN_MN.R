### fit models for simulated data
check.AN = function(M = 10, # Simulation size 
                                n = 50,
                                # Total # of time/new individuals: n = 20                                                #
                                k_0 = 10,
                                # Total # of clusters at baseline (time 0): k_0 = 20
                                Beta00 = 0.5,
                                Beta10 = 0.5,
                                
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
                      # Total # of time/new individuals: n = 20                                                #
                      k_0 = k_0,
                      # Total # of clusters at baseline (time 0): k_0 = 20
                      Beta00 = Beta00,
                      Beta10 = Beta10) 
   
    x_jt_minus_x_new = SN1$x_jt_minus_x_new # dim: [n-(k_t-k_0)] * (k_t) or (k_t-1)
    Z_vec = SN1$Z_vec
    
    XX = plyr::ldply(x_jt_minus_x_new, rbind) ## Cov
    XX = as.data.frame(XX)
    XX = as.matrix(XX)
    
    ############## BN ###############
    
    ZZ.BN = c()
    t_remove = c()
    for (t in 1:length(Z_vec)){
      if(Z_vec[[t]][1] != 1){
        ZZ.BN = rbind(ZZ.BN, c(Z_vec[[t]][-1]))
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
    Beta.est.BN = matrix(0, nrow=k_0, ncol=2)
    for (k in (1:k_0)) {
      MLE2 = glm(c(ZZ.BN[, k]) ~ c(XX[, k]), family = "binomial")
      Beta.est.BN[k, ] = MLE2$coef
      
      l.beta.est.BN = Beta.est.BN[k,1] + Beta.est.BN[k,2] * XX[, k]
      Rho_k.est = exp(l.beta.est.BN)/(1+exp(l.beta.est.BN)) # BERNOULLI RHO
      Rho_k.est.mean[k] = mean(Rho_k.est)
      
    }
    
    # BN for ALL:
    MLE2.all = glm(c(ZZ.BN) ~ c(XX), family = "binomial")
    Beta.est.BN.all = MLE2.all$coef
    Rho.est.BN = as.numeric()
    for(k in (1:k_0)){
      l.beta.est.BN.all = Beta.est.BN.all[1] + Beta.est.BN.all[2] * XX[, k]
      Rho.est.BN[k] = mean(exp(l.beta.est.BN.all)/(1+exp(l.beta.est.BN.all))) # BERNOULLI RHO
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
          ZZ.MN = rbind(ZZ.MN, c(Z_vec[[t]]))
    }
    
    ### Inf loop of MN fit ###
    beta = c(0.5, 0.5)
    MLE1.MN = optim(beta, logit.nll, XX=XX, ZZ=ZZ.MN,
                 method = "BFGS", hessian=TRUE)
    Beta.est.MN = MLE1.MN$par
    Var.est.MN = solve(MLE1.MN$hessian)
    
    ### MN proportions ###
    l.beta.est.MN = Beta.est.MN[1] + Beta.est.MN[2] * XX
    exp.l.beta.MN = exp(l.beta.est.MN)
    p_k.est =  matrix(as.numeric(0), nrow=nrow(XX), ncol=-1+ncol(XX) + 1)
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
