network.sim <-
  function(n = 20,
           # Total # of time/new individuals: n = 20                                                #
           k_0 = 5,
           # Total # of clusters at baseline (time 0): k_0 = 5
           Beta00 = 0.8,
           Beta10 = 0.5) {
    
    n = n
    k_0 = k_0
    Beta00 = Beta00 
    Beta10 = Beta10 
    
    #####################################################################
    s_j0 = rep(1, k_0) # Cluster size at baseline: s_j0 = 1, j=1,...,k_0
    x_j0 = rnorm(k_0, 3*(1:k_0), 0.1)  # explanatory var for the jth cluster at baseline (time 0) (j=1,...,k_0):
    #####################################################################
    
    # To track and update at each t:
    
    ## T=0
    k_t = k_0 
    s_jt = s_j0
    x_jt = x_j0

    x_jt_minus_x_new = Z_vec = list()
    
    for (t in 1:n){
      ind = sample(1:k_0, 1)
      x_new = rnorm(1, 3*ind, 0.01)
      
      ### V2: Measure the distance of NLC with the group center
      xx = sapply(sapply(1/abs(x_jt - x_new), max, 0.001), min, 100)
      l.beta = Beta00 * 1 + Beta10*xx
      upper = c(1, exp(l.beta))
      lower = sum(upper)
      eta_vec = upper/lower
      genMN = rmultinom(1, size = 1, prob = eta_vec)
      j = which.max(genMN)
      
      # To track
      x_jt_minus_x_new = c(x_jt_minus_x_new, list(xx))
      Z_vec = c(Z_vec, list(genMN))
      
      if(j != 1){
        # Parameter Updates:
        #####################################################################   
        k_t.plus.1 = k_t
        s_jt.plus.1 = s_jt
        s_jt.plus.1[j-1] = s_jt[j-1] + 1
        
        x_t.plus.1 = x_jt
        x_t.plus.1[j-1] = mean(c(x_new, x_jt[j-1]))
        
        s_jt = s_jt.plus.1
        k_t = k_t.plus.1
        x_jt = x_t.plus.1
      }
    }
    
    
    return(
      list(
        s_jt = s_jt,
        k_0 = k_0,
        k_t = k_t,
        x_jt_minus_x_new = x_jt_minus_x_new,
        Z_vec = Z_vec
      )
    )
  }

