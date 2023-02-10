rm(list=ls())
source('simulation_functions/perform_jointGHS_simulation_gembag.R')
source('simulation_functions/help_functions.R')


nCores=5
# Perform simulation study

# Choose which scenarios to simulate
run.K2 = F
run.case2 = T
run.case3 = T

# 100 simulations per case
N = 100

# K = 2 data sets, of various similarity -------------
p = 50
K=2
n.vals = c(50,80)
fracs.disagreement = c(0,0.2,0.4,0.6,0.8,1)

if(run.K2){
  # Case 1: datasets from same distribution
  set.seed(1234)
  res.1 = perform_jointGHS_simulation_gembag(K,n.vals, p, N, frac.disagreement = fracs.disagreement[1], nCores = nCores, save.edgeagreement=TRUE)
  
  # Case 2: datasets from similar distributions (80% edge agreement)
  set.seed(1234)
  res.2 = perform_jointGHS_simulation_gembag(K,n.vals, p, N, frac.disagreement = fracs.disagreement[2], nCores = nCores, save.edgeagreement=TRUE)
  
  # Case 3: datasets from slightly related distributions (60% edge agreement)
  set.seed(1234)
  res.3 = perform_jointGHS_simulation_gembag(K,n.vals, p, N, frac.disagreement = fracs.disagreement[3], nCores = nCores, save.edgeagreement=TRUE)
  
  # Case 4: datasets from slightly related distributions (40% edge agreement)
  set.seed(1234)
  res.4 = perform_jointGHS_simulation_gembag(K,n.vals, p, N, frac.disagreement = fracs.disagreement[4], nCores = nCores, save.edgeagreement=TRUE)
  
  # Case 5: datasets from slightly related distributions (20% edge agreement)
  set.seed(1234)
  res.5 = perform_jointGHS_simulation_gembag(K,n.vals, p, N, frac.disagreement = fracs.disagreement[5], nCores = nCores, save.edgeagreement=TRUE)
  
  # Case 6: datasets from unrelated distributions
  set.seed(1234)
  res.6 = perform_jointGHS_simulation_gembag(K,n.vals, p, N, frac.disagreement = fracs.disagreement[6], nCores = nCores, save.edgeagreement=TRUE)
  
  # Save results
  res.K2.gembag = list(res.1, res.2, res.3, res.4, res.5, res.6)
  save(res.K2.gembag, file="data/jointGHS_simulations_K2_GemBag.Rdata")
  
  res.K2.gembag[[1]]$true.sparsity
  
  # Print results
  print_results_jointGHS_show_SD_gembag(res.K2.gembag, fracs.disagreement, show.edgedisagreement = T)
  # Look at edge disagreement (1- fraction of predicted edges that were in common)
}

# K = 2 data sets, of various similarity -------------
p = 100
K=2
n.vals = c(100,150)
fracs.disagreement = c(0,0.2,0.4,0.6,0.8,1)

if(run.case2){
  # Case 1: datasets from same distribution
  set.seed(1234)
  res.1 = perform_jointGHS_simulation_gembag(K,n.vals, p, N, frac.disagreement = fracs.disagreement[1], nCores = nCores, save.edgeagreement=TRUE)
  
  # Case 2: datasets from similar distributions (80% edge agreement)
  set.seed(1234)
  res.2 = perform_jointGHS_simulation_gembag(K,n.vals, p, N, frac.disagreement = fracs.disagreement[2], nCores = nCores, save.edgeagreement=TRUE)
  
  # Case 3: datasets from slightly related distributions (60% edge agreement)
  set.seed(1234)
  res.3 = perform_jointGHS_simulation_gembag(K,n.vals, p, N, frac.disagreement = fracs.disagreement[3], nCores = nCores, save.edgeagreement=TRUE)
  
  # Case 4: datasets from slightly related distributions (40% edge agreement)
  set.seed(1234)
  res.4 = perform_jointGHS_simulation_gembag(K,n.vals, p, N, frac.disagreement = fracs.disagreement[4], nCores = nCores, save.edgeagreement=TRUE)
  
  # Case 5: datasets from slightly related distributions (20% edge agreement)
  set.seed(1234)
  res.5 = perform_jointGHS_simulation_gembag(K,n.vals, p, N, frac.disagreement = fracs.disagreement[5], nCores = nCores, save.edgeagreement=TRUE)
  
  # Case 6: datasets from unrelated distributions
  set.seed(1234)
  res.6 = perform_jointGHS_simulation_gembag(K,n.vals, p, N, frac.disagreement = fracs.disagreement[6], nCores = nCores, save.edgeagreement=TRUE)
  
  # Save results
  res.case2.gembag = list(res.1, res.2, res.3, res.4, res.5, res.6)
  save(res.case2.gembag, file="extended_simulation_studies/data/jointGHS_simulations_case2_GemBag.Rdata")
  
  res.case2.gembag[[1]]$true.sparsity
  
  # Print results
  print_results_jointGHS_show_SD_gembag(res.case2.gembag, fracs.disagreement, show.edgedisagreement = T)
  # Look at edge disagreement (1- fraction of predicted edges that were in common)
}

# K = 2 data sets, of various similarity -------------
p = 100
K=3
n.vals = c(100,150,120)
fracs.disagreement = c(0,0.2,0.4,0.6,0.8,1)

if(run.case3){
  # Case 1: datasets from same distribution
  set.seed(1234)
  res.1 = perform_jointGHS_simulation_gembag(K,n.vals, p, N, frac.disagreement = fracs.disagreement[1], nCores = nCores, save.edgeagreement=TRUE)
  
  # Case 2: datasets from similar distributions (80% edge agreement)
  set.seed(1234)
  res.2 = perform_jointGHS_simulation_gembag(K,n.vals, p, N, frac.disagreement = fracs.disagreement[2], nCores = nCores, save.edgeagreement=TRUE)
  
  # Case 3: datasets from slightly related distributions (60% edge agreement)
  set.seed(1234)
  res.3 = perform_jointGHS_simulation_gembag(K,n.vals, p, N, frac.disagreement = fracs.disagreement[3], nCores = nCores, save.edgeagreement=TRUE)
  
  # Case 4: datasets from slightly related distributions (40% edge agreement)
  set.seed(1234)
  res.4 = perform_jointGHS_simulation_gembag(K,n.vals, p, N, frac.disagreement = fracs.disagreement[4], nCores = nCores, save.edgeagreement=TRUE)
  
  # Case 5: datasets from slightly related distributions (20% edge agreement)
  set.seed(1234)
  res.5 = perform_jointGHS_simulation_gembag(K,n.vals, p, N, frac.disagreement = fracs.disagreement[5], nCores = nCores, save.edgeagreement=TRUE)
  
  # Case 6: datasets from unrelated distributions
  set.seed(1234)
  res.6 = perform_jointGHS_simulation_gembag(K,n.vals, p, N, frac.disagreement = fracs.disagreement[6], nCores = nCores, save.edgeagreement=TRUE)
  
  # Save results
  res.case3.gembag = list(res.1, res.2, res.3, res.4, res.5, res.6)
  save(res.case3.gembag, file="GemBag/data/jointGHS_simulations_case3_GemBag.Rdata")
  
  res.case3.gembag[[1]]$true.sparsity
  
  # Print results
  print_results_jointGHS_show_SD_gembag(res.case3.gembag, fracs.disagreement, show.edgedisagreement = T)
  # Look at edge disagreement (1- fraction of predicted edges that were in common)
}
