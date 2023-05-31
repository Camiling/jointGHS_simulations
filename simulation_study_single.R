source('simulation_functions/perform_fastGHS_simulation.R')
source('simulation_functions/help_functions.R')
source('simulation_functions/perform_time_simulation.R')

# Write print function

nCores = 20 # If using HPC

# Perform simulation study

perform_GHS_sims = FALSE
perform_GHS_sims_redo3 = FALSE # Redo case 3
perform_larger_sims = FALSE
perform_GHS_sims_additional = TRUE

# 20 simulations per case
N = 20

if(perform_GHS_sims){
  # Case 1: p=50,n=100
  set.seed(12345)
  p = 50
  n = 100
  res.1 = perform_fastGHS_simulation(n, p, N, nCores = nCores)

  # Case 2: p=50, n=200
  set.seed(1234)
  p = 50
  n = 200
  res.2 = perform_fastGHS_simulation(n, p, N, nCores = nCores)

  # Case 3: p=100, n=100
  set.seed(1234)
  p = 100
  n = 100
  res.3 = perform_fastGHS_simulation(n, p, N, nCores = nCores)

  # Case 4: p=100, n=200
  set.seed(1234)
  p = 100
  n = 200
  res.4 = perform_fastGHS_simulation(n, p, N, nCores = nCores)
  
  # Save results
  res.fast = list(res.1, res.2, res.3, res.4)
  #save(res.fast, file="data/fastGHS_simulations.Rdata")
  load("data/fastGHS_simulations.Rdata")
  
  # Print results
  print_results_fastGHS(res.fast, show.interval=F, show.sd=T)
}

if(perform_larger_sims){ # Not used
  # Case 5: p=200,n=150 (Here GHS with Gibbs fails, so we do not run)
  set.seed(12345)
  p = 200
  n = 150
  res.5 = perform_fastGHS_simulation(n, p, N, include.GHS = F, nCores = nCores)
 
  # Case 6: p=400, n=300
  set.seed(1234)
  p = 300
  n = 300
  res.6 = perform_fastGHS_simulation(n, p, N, include.GHS = F, nCores = nCores)
  
  # Save results
  res.fast = list(res.5,res.6)
  #save(res.fast, file="data/fastGHS_simulations_large.Rdata")
  load("data/fastGHS_simulations_large.Rdata")
  
  print_results_fastGHS(res.fast, show.interval=F, show.sd=T, include.GHS = F)
}

if(perform_GHS_sims_redo3){
  # Case 3: p=100, n=100
  set.seed(12345)
  p = 100
  n = 150
  res.3 = perform_fastGHS_simulation(n, p, N, nCores = nCores)
  #save(res.3, file="data/fastGHS_simulations_redo3.Rdata")
  load("data/fastGHS_simulations.Rdata")
  load("data/fastGHS_simulations_redo3.Rdata")
  res.fast[[3]] = res.3
  # Print results
  print_results_fastGHS(res.fast, show.interval=F, show.sd=T)
}

if(perform_GHS_sims_additional){
  # A1: p=50, n=50
  set.seed(1234)
  p = 50
  n = 50
  res.a1 = perform_fastGHS_simulation(n, p, N, nCores = nCores)
  
  # A2: p=50, n=150
  set.seed(1234)
  p = 50
  n = 150
  res.a2 = perform_fastGHS_simulation(n, p, N, nCores = nCores)
  
  # A3: p=100, n=250
  set.seed(1234)
  p = 100
  n = 250
  res.a3 = perform_fastGHS_simulation(n, p, N, nCores = nCores)
  
  # Save results
  res.a = list(res.a1, res.a2, res.a3)
  #save(res.a, file="data/fastGHS_simulations_additional.Rdata")
  load("data/fastGHS_simulations_additional.Rdata")
  print_results_fastGHS(res.a, show.interval=F, show.sd=T)
}



  