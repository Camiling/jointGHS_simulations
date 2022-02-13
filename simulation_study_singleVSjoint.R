source('simulation_functions/perform_jointGHS_simulation.R')
source('simulation_functions/help_functions.R')

nCores = 40 # If using HPC

# Perform simulation study comparing single-network fastGHS to jointGHS
# We force fastGHS to the same sparsity as jointGHS for better comparison

# Choose which scenarios to simulate
run.K2 = T
run.K4 = F
run.K10 = F

# 20 simulations per case
N = 40 # CHANGED TO 40 for new7, with new5 seed and 80, 100 n vals
p = 50

# K = 2 data sets, of various similarity -------------
K=2
seed.1 = 12345 
n.vals=c(80,100) 
fracs.disagreement = c(0,0.2,0.4,0.6,0.8,1)

if(run.K2){
  # Case 1: datasets from same distribution
  set.seed(seed.1)
  res.1 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = fracs.disagreement[1],include.SSJGL=FALSE, include.JGL = FALSE, nCores = nCores, singleVSjoint=T)
  
  # Case 2: datasets from similar distributions (80% edge agreement)
  set.seed(seed.1)
  res.2 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = fracs.disagreement[2], include.SSJGL=FALSE, include.JGL = FALSE,nCores = nCores, singleVSjoint=T)
  
  # Case 3: datasets from slightly related distributions (60% edge agreement)
  set.seed(seed.1)
  res.3 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = fracs.disagreement[3],include.SSJGL=FALSE, include.JGL = FALSE, nCores = nCores, singleVSjoint=T)
  
  # Case 4: datasets from slightly related distributions (40% edge agreement)
  set.seed(seed.1)
  res.4 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = fracs.disagreement[4], include.SSJGL=FALSE, include.JGL = FALSE,nCores = nCores, singleVSjoint=T)
  
  # Case 5: datasets from slightly related distributions (20% edge agreement)
  set.seed(seed.1)
  res.5 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = fracs.disagreement[5], include.SSJGL=FALSE, include.JGL = FALSE, nCores = nCores, singleVSjoint=T)
  
  # Case 6: datasets from unrelated distributions
  set.seed(seed.1)
  res.6 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = fracs.disagreement[6], include.SSJGL=FALSE, include.JGL = FALSE, nCores = nCores, singleVSjoint=T)
  
  # Save results
  res.K2 = list(res.1, res.2, res.3, res.4, res.5, res.6)
  save(res.K2, file="data/jointGHS_simulations_K2_singleVSjoint_final4.Rdata")
  
  # Print results
  print_results_jointGHS(res.K2, fracs.disagreement, show.interval=F, show.sd=F, include.SSJGL = F, include.JGL = F)
}


# K = 4 data sets, of various similarity -------------
K.2 = 4
n.vals.2 = c(80, 100, 100, 100) 

if(run.K4){
  # Case 7: datasets from same distribution
  set.seed(12345)
  res.7 = perform_jointGHS_simulation(K.2,n.vals.2, p, N, frac.disagreement = fracs.disagreement[1], include.SSJGL=FALSE, include.JGL = FALSE,nCores = nCores, singleVSjoint=T)
  
  # Case 8: datasets from similar distributions (80% edge agreement)
  set.seed(12345)
  res.8 = perform_jointGHS_simulation(K.2,n.vals.2, p, N, frac.disagreement = fracs.disagreement[2],include.SSJGL=FALSE, include.JGL = FALSE, nCores = nCores, singleVSjoint=T)
  
  # Case 8: datasets from slightly related distributions (60% edge agreement)
  set.seed(12345)
  res.9 = perform_jointGHS_simulation(K.2,n.vals.2, p, N, frac.disagreement = fracs.disagreement[3],include.SSJGL=FALSE, include.JGL = FALSE, nCores = nCores, singleVSjoint=T)
  
  # Case 10: datasets from slightly related distributions (40% edge agreement)
  set.seed(12345)
  res.10 = perform_jointGHS_simulation(K.2,n.vals.2, p, N, frac.disagreement = fracs.disagreement[4],include.SSJGL=FALSE, include.JGL = FALSE, nCores = nCores, singleVSjoint=T)
  
  # Case 11: datasets from slightly related distributions (20% edge agreement)
  set.seed(12345)
  res.11 = perform_jointGHS_simulation(K.2,n.vals.2, p, N, frac.disagreement = fracs.disagreement[5],include.SSJGL=FALSE, include.JGL = FALSE, nCores = nCores, singleVSjoint=T)
  
  # Case 12: datasets from unrelated distributions
  set.seed(12345)
  res.12 = perform_jointGHS_simulation(K.2,n.vals.2, p, N, frac.disagreement = fracs.disagreement[6],include.SSJGL=FALSE, include.JGL = FALSE, nCores = nCores, singleVSjoint=T)
  
  # Save results
  res.K4.list = list(res.7, res.8, res.9, res.10, res.11, res.12)
  save(res.K4.list, file="data/jointGHS_simulations_K4_singleVSjoint_final.Rdata")
  
  # Print results
  print_results_jointGHS(res.K4.list,fracs.disagreement, show.interval=F,show.sd=F,include.SSJGL = F, include.JGL = F)
}



# K = 10 data sets, of various similarity -------------
K.3=10
n.vals.3 = c(rep(80,5), rep(100,5))
seed.10 = 1234 


if(run.K10){
  # Case 13: datasets from same distribution
  set.seed(seed.10)
  res.13 = perform_jointGHS_simulation(K.3,n.vals.3, p, N, frac.disagreement = fracs.disagreement[1],include.SSJGL = F, include.JGL = F, nCores = nCores, singleVSjoint=T)
  
  # Case 14: datasets from similar distributions (80% edge agreement)
  set.seed(seed.10)
  res.14 = perform_jointGHS_simulation(K.3,n.vals.3, p, N, frac.disagreement = fracs.disagreement[2],include.SSJGL = F, include.JGL = F, nCores = nCores, singleVSjoint=T)
  
  # Case 15: datasets from slightly related distributions (60% edge agreement)
  set.seed(seed.10)
  res.15 = perform_jointGHS_simulation(K.3,n.vals.3, p, N, frac.disagreement = fracs.disagreement[3],include.SSJGL = F, include.JGL = F, nCores = nCores, singleVSjoint=T)
  
  # Case 16: datasets from slightly related distributions (40% edge agreement)
  set.seed(seed.10)
  res.16 = perform_jointGHS_simulation(K.3,n.vals.3, p, N, frac.disagreement = fracs.disagreement[4],include.SSJGL = F, include.JGL = F,nCores = nCores, singleVSjoint=T)
  
  # Case 17: datasets from slightly related distributions (20% edge agreement)
  set.seed(seed.10)
  res.17 = perform_jointGHS_simulation(K.3,n.vals.3, p, N, frac.disagreement = fracs.disagreement[5], include.SSJGL = F, include.JGL = F,nCores = nCores, singleVSjoint=T)
  
  # Case 18: datasets from unrelated distributions
  set.seed(seed.10)
  res.18 = perform_jointGHS_simulation(K.3,n.vals.3, p, N, frac.disagreement = fracs.disagreement[6], include.SSJGL = F, include.JGL = F,nCores = nCores, singleVSjoint=T)
  
  # Save results
  res.K10.list = list(res.13, res.14, res.15, res.16, res.17, res.18)  
  save(res.K10.list, file="data/jointGHS_simulations_K10_singleVSjoint_final.Rdata")
  
  # Print results
  print_results_jointGHS(res.K10.list,fracs.disagreement, show.interval=F, show.sd=T, include.SSJGL = F, include.JGL = F,collapse=T)
  
}



