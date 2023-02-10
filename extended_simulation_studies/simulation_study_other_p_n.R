rm(list=ls())
source('simulation_functions/perform_jointGHS_simulation.R')
source('simulation_functions/help_functions.R')

nCores = 56 # If using HPC

# Perform simulation study

# Choose which scenarios to simulate
run.case1 = F
run.case2 = T

# 100 simulations per case
N = 100

# K = 2 data sets, p=100 of various similarity -------------
p=100
K=2
n.vals = c(100,150)
fracs.disagreement = c(0,0.2,0.4,0.6,0.8,1)

if(run.case1){
  # Case 1: datasets from same distribution
  set.seed(1234)
  res.1 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = fracs.disagreement[1], nCores = nCores, save.edgeagreement=T)
  save(res.1, file="extended_simulation_studies/data/jointGHS_simulations_case1_1.Rdata")
  
  # Case 2: datasets from similar distributions (80% edge agreement)
  set.seed(1234)
  res.2 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = fracs.disagreement[2], nCores = nCores, save.edgeagreement=T)
  save(res.2, file="extended_simulation_studies/data/jointGHS_simulations_case1_2.Rdata")
  
  # Case 3: datasets from slightly related distributions (60% edge agreement)
  set.seed(1234)
  res.3 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = fracs.disagreement[3], nCores = nCores, save.edgeagreement=T)
  save(res.3, file="extended_simulation_studies/data/jointGHS_simulations_case1_3.Rdata")
  
  # Case 4: datasets from slightly related distributions (40% edge agreement)
  set.seed(1234)
  res.4 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = fracs.disagreement[4], nCores = nCores, save.edgeagreement=T)
  save(res.4, file="extended_simulation_studies/data/jointGHS_simulations_case1_4.Rdata")
  
  # Case 5: datasets from slightly related distributions (20% edge agreement)
  set.seed(1234)
  res.5 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = fracs.disagreement[5], nCores = nCores, save.edgeagreement=T)
  save(res.5, file="extended_simulation_studies/data/jointGHS_simulations_case1_5.Rdata")
  
  # Case 6: datasets from unrelated distributions
  set.seed(1234)
  res.6 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = fracs.disagreement[6], nCores = nCores, save.edgeagreement=T)
  save(res.6, file="extended_simulation_studies/data/jointGHS_simulations_case1_6.Rdata")
  
  # Save results
  res.case1 = list(res.1, res.2, res.3, res.4, res.5, res.6)
  save(res.case1, file="extended_simulation_studies/data/jointGHS_simulations_case1.Rdata")
  
  # Print results
  #print_results_jointGHS_show_SD(res.case1, fracs.disagreement,include.GHS = F, show.edgedisagreement = T)
}

# K = 3 data sets, p=100 of various similarity -------------
p=100
K=3
n.vals = c(100,150,120)
fracs.disagreement = c(0,0.2,0.4,0.6,0.8,1)

if(run.case2){
  # Case 1: datasets from same distribution
  set.seed(1234)
  res.1 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = fracs.disagreement[1], nCores = nCores, save.edgeagreement=T)
  save(res.1, file="extended_simulation_studies/data/jointGHS_simulations_case2_1.Rdata")
  
  # Case 2: datasets from similar distributions (80% edge agreement)
  set.seed(1234)
  res.2 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = fracs.disagreement[2], nCores = nCores, save.edgeagreement=T)
  save(res.2, file="extended_simulation_studies/data/jointGHS_simulations_case2_2.Rdata")
  
  # Case 3: datasets from slightly related distributions (60% edge agreement)
  set.seed(1234)
  res.3 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = fracs.disagreement[3], nCores = nCores, save.edgeagreement=T)
  save(res.3, file="extended_simulation_studies/data/jointGHS_simulations_case2_3.Rdata")
  
  # Case 4: datasets from slightly related distributions (40% edge agreement)
  set.seed(1234)
  res.4 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = fracs.disagreement[4], nCores = nCores, save.edgeagreement=T)
  save(res.4, file="extended_simulation_studies/data/jointGHS_simulations_case2_4.Rdata")
  
  # Case 5: datasets from slightly related distributions (20% edge agreement)
  set.seed(1234)
  res.5 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = fracs.disagreement[5], nCores = nCores, save.edgeagreement=T)
  save(res.5, file="extended_simulation_studies/data/jointGHS_simulations_case2_5.Rdata")
  
  # Case 6: datasets from unrelated distributions
  set.seed(1234)
  res.6 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = fracs.disagreement[6], nCores = nCores, save.edgeagreement=T)
  save(res.6, file="extended_simulation_studies/data/jointGHS_simulations_case2_6.Rdata")
  
  # Save results
  res.case2 = list(res.1, res.2, res.3, res.4, res.5, res.6)
  save(res.case2, file="extended_simulation_studies/data/jointGHS_simulations_case2.Rdata")
  
  # Print results
  #print_results_jointGHS_show_SD(res.case2, fracs.disagreement,include.GHS = F, show.edgedisagreement = T)
}


