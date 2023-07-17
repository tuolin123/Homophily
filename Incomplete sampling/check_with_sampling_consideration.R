
### Check nomianal level of alpha
check.AN = function(M = 10, # Simulation size 
                    n = 50,
                    # Total # of time/new individuals: n = 20                                                #
                    k_0 = 5,
                    # Total # of clusters at baseline (time 0): k_0 = 20
                    k_tot = 10,
                    Beta00 = 0.5,
                    Beta10 = 0.5,
                    Beta20 = 0.5,
                    size1 = 300,
                    size2 = 1100,
                    midsize = 200) {
  
  Beta0_chisq.MN = Beta1_chisq.MN  = c()
  Beta_est.BN = Var_est.BN = list()
  Beta_est.MN = Var_est.MN = list()
  
  for (m in 1:M) {
    
    ### Generating Network ###
    SN1 = network.sim(n = n,
                      # Total # of time/new individuals: n = 20                                                #
                      k_0 = k_0,
                      # Total # of clusters at baseline (time 0): k_0 = 20
                      k_tot = k_tot,
                      Beta00 = Beta00,
                      Beta10 = Beta10,
                      Beta20 = Beta20,
                      size1 = size1,
                      size2 = size2,
                      midsize = midsize) 
    ind = which(SN1$x_his == 0)
    ind_exclude = sample(ind, floor(0.5*length(ind)))
    ind_sample = (1:length(SN1$x_prop))[-ind_exclude]
    
    x1_jt_minus_x_new = SN1$x1_jt_minus_x_new[ind_sample]
    x2_jt_minus_x_new = SN1$x2_jt_minus_x_new[ind_sample]
    Z_vec = SN1$Z_vec[ind_sample]
    
    XX1 = plyr::ldply(x1_jt_minus_x_new, rbind) ## Cov
    XX1 = as.data.frame(XX1)
    XX1 = as.matrix(XX1)
    
    XX2 = plyr::ldply(x2_jt_minus_x_new, rbind) ## Cov
    XX2 = as.data.frame(XX2)
    XX2 = as.matrix(XX2)
    
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
      XX1 = XX1
      XX2 = XX2
    }else{
      XX1 = XX1[-t_remove, ]
      XX2 = XX2[-t_remove, ]
    }
    
    ### Inf loop for BN
    # MLE1.BN = glm(c(ZZ.BN) ~ c(XX), family = "binomial")
    # Beta.est.BN = MLE1.BN$coef
    # Var.est.BN = diag((summary(MLE1.BN)$coefficients[,2])^2)
    Rho_k.est.mean =  as.numeric()
    Beta.est.BN = matrix(0, nrow=ncol(XX1), ncol=3)
    for (k in (1:ncol(XX1))) {
      datfit = data.frame(cbind(ZZ.BN[, k], XX1[, k], XX2[, k]))
      datfit = datfit[complete.cases(datfit),]
      colnames(datfit) = c("y", "x1", "x2")
      MLE2 = glm(y ~ x1 + x2, data = datfit, family = "binomial")
      Beta.est.BN[k, ] = MLE2$coef
      
      l.beta.est.BN = Beta.est.BN[k,1] + Beta.est.BN[k,2] * XX1[, k] + 
        Beta.est.BN[k,3] * XX2[, k]
      Rho_k.est = exp(l.beta.est.BN)/(1+exp(l.beta.est.BN)) # BERNOULLI RHO
      Rho_k.est.mean[k] = mean(Rho_k.est, na.rm = T)
      
    }
    
    #	Check if the p can be reversed back to be the same:
    
    # BN for ALL:
    datfit = data.frame(cbind(c(ZZ.BN), c(XX1), c(XX2)))
    datfit = datfit[complete.cases(datfit),]
    colnames(datfit) = c("y", "x1", "x2")
    
    MLE2.all = glm(y ~ x1 + x2, data = datfit, family = "binomial")
    Beta.est.BN.all = MLE2.all$coef
    Rho.est.BN = as.numeric()
    for(k in (1:ncol(XX1))){
      l.beta.est.BN.all = Beta.est.BN.all[1] + Beta.est.BN.all[2] * XX1[, k] +
        Beta.est.BN.all[3] * XX2[, k]
      Rho.est.BN[k] = mean(exp(l.beta.est.BN.all)/(1+exp(l.beta.est.BN.all)), na.rm = T) # BERNOULLI RHO
    }
    
    #BN.p_k.est = round(colMeans(Rho.est.BN, na.rm = T),4)
    BN.p_k.est = round(Rho.est.BN, 4)
    
    Beta_est.BN[[m]] = as.vector(Beta.est.BN.all)
    #Var_est.BN[[m]] = (summary(MLE1.BN)$coefficients[,2])^2
    
    
    ############## MN ###############
    XX1 = plyr::ldply(x1_jt_minus_x_new, rbind) ## Cov
    XX1 = as.data.frame(XX1)
    XX1 = as.matrix(XX1)
    
    XX2 = plyr::ldply(x2_jt_minus_x_new, rbind) ## Cov
    XX2 = as.data.frame(XX2)
    XX2 = as.matrix(XX2)
    
    ZZ.MN = c()
    for (t in 1:length(Z_vec)){
      update = as.vector(Z_vec[[t]])
      length(update) = max(sapply(Z_vec, length))
      ZZ.MN = rbind(ZZ.MN, update)
    }
    ZZ.MN = as.matrix(ZZ.MN)
    
    ### Inf loop of MN fit ###
    beta = c(0.5, 0.5, 0.5)
    vec_len = sapply(Z_vec, length)
    MLE1.MN = optim(beta, logit.nll.all, XX1=XX1, XX2=XX2, ZZ=ZZ.MN,
                    len = vec_len, ind = 6:11,
                    method = "BFGS", hessian=TRUE)
    
    
    Beta.est.MN = MLE1.MN$par
    #Var.est.MN = solve(MLE1.MN$hessian)
    
    ### MN proportions ###
    l.beta.est.MN = Beta.est.MN[1] + Beta.est.MN[2] * XX1 + 
      Beta.est.MN[3] * XX2
    exp.l.beta.MN = exp(l.beta.est.MN)
    p_k.est =  matrix(as.numeric(0), nrow=nrow(XX1), ncol=ncol(XX1))
    for(i in (1:nrow(XX1))){
      col_k = length(exp.l.beta.MN[i,])
      #p_k.est[i,] = c(exp.l.beta[i,],1)/(1+sum(exp.l.beta[i,],na.rm = T))
      #p_k.est[i,] = exp.l.beta.MN[i,]/(1+sum(exp.l.beta.MN[i,],na.rm = T))
      p_k.est[i,] = exp.l.beta.MN[i,]/(1+sum(exp.l.beta.MN[i,],na.rm = T))
      
    }
    MN.p_k.est = round(colMeans(p_k.est, na.rm=T),4)
    
    
    #Beta0_chisq.MN = c(Beta0_chisq.MN, chisq_stat(Beta0_matrix, (Beta.est.MN-Beta00), Var.est.MN))
    #Beta1_chisq.MN = c(Beta1_chisq.MN, chisq_stat(Beta1_matrix, (Beta.est.MN-Beta10), Var.est.MN))
    
    Beta_est.MN[[m]] = MLE1.MN$par
    #Var_est.MN[[m]] = diag(solve(MLE1.MN$hessian))
    
    
  }
  
  Beta.emp.BN = rowMeans(as.data.frame(Beta_est.BN))
  Beta.emp.MN = rowMeans(as.data.frame(Beta_est.MN))
  sd.emp.MN = apply(as.data.frame(Beta_est.MN), 1, sd)
  
  
  return(list( 
    
    ### BN ###
    Rho_k.est.mean = Rho_k.est.mean, # BN for each column
    BN.p_k.est = BN.p_k.est,         # BN for ALL
    Beta.emp.BN = Beta.emp.BN,
    
    
    ### MN ###
    MN.p_k.est = MN.p_k.est,
    Beta.emp.MN = Beta.emp.MN,
    sd.emp.MN = sd.emp.MN
  ))
}
