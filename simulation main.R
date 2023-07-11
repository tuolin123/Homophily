## Scenario I: Simulation study
# Fixed cluster number
source('Scenario I/Non_increasing_cluster_number/Sim_network_Func_diffX.R')
source('Scenario I/Non_increasing_cluster_number/Optim_func.R')
source('Scenario I/Non_increasing_cluster_number/check_AN_BN_MN.R')

Beta0_matrix = matrix(c(1, 0), nrow = 1)
Beta1_matrix = matrix(c(0, 1), nrow = 1)

n= 1000

T0 = check.AN(M = 500, # Simulation size
              n = n,
              # Total # of time/new individuals: n = 1000                                                #
              k_0 = 10,
              # Total # of clusters at baseline (time 0): k_0 = 10
              Beta00 = -0.8,
              Beta10 = 0.1,
              
              Beta0_matrix,
              Beta1_matrix)

# Increasing cluster number
source('Scenario I/Increasing_cluster_number/Sim_network_Func_Increase.R')
source('Scenario I/Increasing_cluster_number/check_AN_BN_MN_increasing.R')
T1 = check.AN(M = 500, # Simulation size 
              n = 1500,
              # Total # of time/new individuals: n = 20                                                #
              k_0 = 5,
              # Total # of clusters at baseline (time 0): k_0 = 5
              k_tot = 10,
              Beta00 = -0.8,
              Beta10 = 0.1,
              size1 = 300,
              size2 = 1100,
              midsize = 200,
              
              Beta0_matrix, 
              Beta1_matrix)

# check the estimates in the independent logistic regression converge
source('Scenario I/Check_independent_logistic_regression/Sim_network_Func_Increase_BN.R')
T2 = check.AN(M = 500, # Simulation size 
              n = 60,
              # Total # of time/new individuals: n = 60                                                #
              k_0 = 5,
              # Total # of clusters at baseline (time 0): k_0 = 5
              k_tot = 10,
              Beta00 = -0.8,
              Beta10 = 0.1,
              size1 = 10,
              size2 = 50,
              midsize = 10,
              
              Beta0_matrix, 
              Beta1_matrix)

T3 = check.AN(M = 500, # Simulation size 
              n = 120,
              # Total # of time/new individuals: n = 20                                                #
              k_0 = 5,
              # Total # of clusters at baseline (time 0): k_0 = 5
              k_tot = 10,
              Beta00 = -0.8,
              Beta10 = 0.1,
              size1 = 20,
              size2 = 100,
              midsize = 20,
              
              Beta0_matrix, 
              Beta1_matrix)


## Scenario II: Real study based simulation
# Fixed cluster number
source('Scenario II/Non_increasing_cluster_number/Sim_network_Func_diffX_real.R')
source('Scenario II/Non_increasing_cluster_number/check_AN_BN_MN_real.R')
source('Scenario II/Non_increasing_cluster_number/Optim_func_real.R')

TR1 = check.AN(M = 500, # Simulation size 
               n = 1000,
               # Total # of time/new individuals: n = 20                                                #
               k_0 = 5,
               # Total # of clusters at baseline (time 0): k_0 = 5
               Beta00 = 1,
               Beta10 = log(0.86),
               Beta20 = log(2.44))

# Increasing cluster number
source('Scenario II/Increasing_cluster_number/Sim_network_Func_Increase_real.R')
source('Scenario II/Increasing_cluster_number/check_AN_BN_MN_increasing_real.R')

TR2 = check.AN(M = 500, # Simulation size 
              n = 1500,
              # Total # of time/new individuals: n = 20                                                #
              k_0 = 5,
              # Total # of clusters at baseline (time 0): k_0 = 5
              k_tot = 10,
              Beta00 = 1,
              Beta10 = log(0.86),
              Beta20 = log(2.44),
              size1 = 300,
              size2 = 1100,
              midsize = 200)

# Check for incomplete sampling case
source('Incomplete sampling/Sim_with_sampling_consideration.R')
source('Incomplete sampling/check_with_sampling_consideration.R')
TRs = check.AN(M = 500, # Simulation size 
              n = 1500,
              # Total # of time/new individuals: n = 20                                                #
              k_0 = 5,
              # Total # of clusters at baseline (time 0): k_0 = 5
              k_tot = 10,
              Beta00 = 1,
              Beta10 = log(0.86),
              Beta20 = log(2.44),
              size1 = 300,
              size2 = 1100,
              midsize = 200)
