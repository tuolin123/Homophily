network.sim <-
  function(n = 20,
           # Total # of time/new individuals: n = 20                                                #
           k_0 = 5,
           # Total # of clusters at baseline (time 0): k_0 = 5
           Beta00 = 0.8,
           Beta10 = 0.5,
           Beta20) {
    
    n = n
    k_0 = k_0
    Beta00 = Beta00 
    Beta10 = Beta10 
    
    #####################################################################
    s_j0 = rep(1, k_0) # Cluster size at baseline: s_j0 = 1, j=1,...,k_0
    meann = c(1968, 1971, 1975, 1978, 1981, 1955, 1961, 1965, 1985, 1991)
    propp = c(0.25, 0.3, 0.35, 0.4, 0.45, 0.1, 0.15, 0.2, 0.5, 0.55)
    sdd = sqrt(1)
    if(k_0 == 5){
      x1_j0 = c(rnorm(1, meann[1], sdd),rnorm(1, meann[2], sdd),rnorm(1, meann[3], sdd),
                      rnorm(1, meann[4], sdd),rnorm(1, meann[5], sdd))
      r2_num_j0 = propp[1:5]*20
      r2_denom_j0 = rep(20,5)
    }
    if(k_0 == 10){
      x1_j0 = rnorm(10, meann, sdd)
      r2_num_j0 = propp[1:10]*20
      r2_denom_j0 = rep(20,10)
    }
    
    #####################################################################
    
    # To track and update at each t:
    
    ## T=0
    k_t = k_0 
    s_jt = s_j0
    x1_jt = x1_j0
    r2_num_jt = r2_num_j0
    r2_denom_jt = r2_denom_j0
    n_clus = rep(1,k_0)
    
    x1_jt_minus_x_new = x2_jt_minus_x_new = Z_vec = list()
    
    for (t in 1:n){
      ind = sample(1:k_0, 1)
      x1_new = round(rnorm(1, meann[ind], sdd),0)
      x2_new = rbinom(1, 1, r2_num_jt[ind]/r2_denom_jt[ind])
      
      ### V2: Measure the distance of NLC with the group center
      xx1 = abs(x1_jt - x1_new)
      r2_jt.plus.1 = (r2_num_jt+1)/(r2_denom_jt+1)
      xx2 = r2_jt.plus.1^x2_new*(1-r2_jt.plus.1)^(1-x2_new)
      l.beta = Beta00 * 1 + Beta10*xx1 + Beta20*xx2
      upper = c(1, exp(l.beta))
      lower = sum(upper)
      eta_vec = upper/lower
      genMN = rmultinom(1, size = 1, prob = eta_vec)
      j = which.max(genMN)
      
      # To track
      x1_jt_minus_x_new = c(x1_jt_minus_x_new, list(xx1))
      x2_jt_minus_x_new = c(x2_jt_minus_x_new, list(xx2))
      Z_vec = c(Z_vec, list(genMN))
      
      ### The first group for multinomial is not group to any of existing group
      if(j != 1){
        # Parameter Updates:
        x1_t.plus.1 = x1_jt
        x1_t.plus.1[j-1] = (x1_new + n_clus[j-1]*x1_jt[j-1])/(1+n_clus[j-1])
        x1_jt = x1_t.plus.1
        
        if(x2_new == 1){
          r2_num_jt[j-1] = r2_num_jt[j-1]+1
        }else{
          r2_num_jt[j-1] = r2_num_jt[j-1]
        }
        r2_denom_jt[j-1] = r2_denom_jt[j-1]+1
        
        n_clus[j-1] = n_clus[j-1] + 1
      } 
    }
    
    
    return(
      list(
        x1_jt_minus_x_new = x1_jt_minus_x_new,
        x2_jt_minus_x_new = x2_jt_minus_x_new,
        Z_vec = Z_vec
      )
    )
  }

