source('simulation_functions/perform_jointGHS_simulation.R')
source('simulation_functions/help_functions.R')

nCores = 56 # If using HPC

# Perform simulation study

# Choose which scenarios to simulate
run.K2 = T

# 100 simulations per case
N = 100
p = 50

# K = 2 data sets, of various similarity -------------
K=2
n.vals = c(50,80)
fracs.disagreement = c(0,0.2,0.4,0.6,0.8,1)

if(run.K2){
  # Case 1: datasets from same distribution
  set.seed(1234)
  res.1 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = fracs.disagreement[1], nCores = nCores)
  cat('ok \n')
  
  # Case 2: datasets from similar distributions (80% edge agreement)
  set.seed(1234)
  res.2 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = fracs.disagreement[2], nCores = nCores)
  
  # Case 3: datasets from slightly related distributions (60% edge agreement)
  set.seed(1234)
  res.3 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = fracs.disagreement[3], nCores = nCores)
  
  # Case 4: datasets from slightly related distributions (40% edge agreement)
  set.seed(1234)
  res.4 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = fracs.disagreement[4], nCores = nCores)
  
  # Case 5: datasets from slightly related distributions (20% edge agreement)
  set.seed(1234)
  res.5 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = fracs.disagreement[5], nCores = nCores)
  
  # Case 6: datasets from unrelated distributions
  set.seed(1234)
  res.6 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = fracs.disagreement[6], nCores = nCores)
  
  # Save results
  res.K2 = list(res.1, res.2, res.3, res.4, res.5, res.6)
  save(res.K2, file="data/jointGHS_simulations_K2.Rdata")
  
  res.K2[[1]]$true.sparsity
  
  # Print results
  print_results_jointGHS(res.K2, fracs.disagreement, show.interval=F, show.sd=T, include.GHS = F)
}




