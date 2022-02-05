source('simulation_functions/perform_jointGHS_simulation.R')
source('simulation_functions/help_functions.R')

nCores = 56 # If using HPC

# Perform simulation study

# Choose which scenarios to simulate
run.K2 = T
run.K4 = T
run.K10 = T

# 100 simulations per case
N = 100
p = 50

# K = 2 data sets, of various similarity -------------
K=2
n.vals = c(50,80)
fracs.disagreement = c(0,0.2,0.4,0.6,1)

if(run.K2){
  # Case 1: datasets from same distribution
  set.seed(1234)
  res.1 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = fracs.disagreement[1], nCores = nCores)
  
  # Case 2: datasets from similar distributions (80% edge agreement)
  set.seed(1234)
  res.2 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = fracs.disagreement[2], nCores = nCores)
  
  # Case 3: datasets from slightly related distributions (60% edge agreement)
  set.seed(1234)
  res.3 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = fracs.disagreement[3], nCores = nCores)
  
  # Case 4: datasets from slightly related distributions (40% edge agreement)
  set.seed(1234)
  res.4 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = fracs.disagreement[4], nCores = nCores)
  
  # Case 5: datasets from unrelated distributions
  set.seed(1234)
  res.5 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = fracs.disagreement[5], nCores = nCores)
  
  # Save results
  res.K2 = list(res.1, res.2, res.3, res.4, res.5)
  save(res.K2, file="data/jointGHS_simulations_K2.Rdata")
  
  # Print results
  #print_results_jointGHS(res.K2, fracs.disagreement, show.interval=F, show.sd=F)
}


# Fewer simulations now
N = 5


# K = 4 data sets, of various similarity -------------
K.2 = 4
n.vals.2 = c(50, 80, 50, 80)


if(run.K4){
  # Case 6: datasets from same distribution
  set.seed(12345)
  res.6 = perform_jointGHS_simulation(K.2,n.vals.2, p, N, frac.disagreement = fracs.disagreement[1], include.GHS=FALSE, include.SSJGL=FALSE, include.JGL = FALSE,nCores = nCores)
  
  # Case 7: datasets from similar distributions (80% edge agreement)
  set.seed(12345)
  res.7 = perform_jointGHS_simulation(K.2,n.vals.2, p, N, frac.disagreement = fracs.disagreement[2],include.GHS=FALSE,include.SSJGL=FALSE, include.JGL = FALSE, nCores = nCores)
  
  # Case 8: datasets from slightly related distributions (60% edge agreement)
  set.seed(12345)
  res.8 = perform_jointGHS_simulation(K.2,n.vals.2, p, N, frac.disagreement = fracs.disagreement[3],include.GHS=FALSE,include.SSJGL=FALSE, include.JGL = FALSE, nCores = nCores)
  
  # Case 9: datasets from slightly related distributions (40% edge agreement)
  set.seed(12345)
  res.9 = perform_jointGHS_simulation(K.2,n.vals.2, p, N, frac.disagreement = fracs.disagreement[4],include.GHS=FALSE,include.SSJGL=FALSE, include.JGL = FALSE, nCores = nCores)
  
  # Case 10: datasets from unrelated distributions
  set.seed(12345)
  res.10 = perform_jointGHS_simulation(K.2,n.vals.2, p, N, frac.disagreement = fracs.disagreement[5],include.GHS=FALSE,include.SSJGL=FALSE, include.JGL = FALSE, nCores = nCores)
  
  
  # Save results
  res.K4.list = list(res.6, res.7, res.8, res.9, res.10)
  save(res.K4.list, file="data/jointGHS_simulations_K4.Rdata")
  
  # Print results
  #print_results_jointGHS(res.K4.list,fracs.disagreement, show.interval=F,include.GHS=FALSE, show.sd=F,include.SSJGL = F, include.JGL = F)
}



# K = 10 data sets, of various similarity -------------
K.3=10
n.vals.3 = c(rep(80,5), rep(100,5))


if(run.K10){
  # Case 11: datasets from same distribution
  set.seed(123)
  res.11 = perform_jointGHS_simulation(K.3,n.vals.3, p, N, frac.disagreement = fracs.disagreement[1],include.GHS=FALSE,include.SSJGL = F, include.JGL = F, nCores = nCores)
  
  # Case 12: datasets from similar distributions (80% edge agreement)
  set.seed(123)
  res.12 = perform_jointGHS_simulation(K.3,n.vals.3, p, N, frac.disagreement = fracs.disagreement[2],include.GHS=FALSE,include.SSJGL = F, include.JGL = F, nCores = nCores)
  
  # Case 13: datasets from slightly related distributions (60% edge agreement)
  set.seed(123)
  res.13 = perform_jointGHS_simulation(K.3,n.vals.3, p, N, frac.disagreement = fracs.disagreement[3],include.GHS=FALSE,include.SSJGL = F, include.JGL = F, nCores = nCores)
  
  # Case 14: datasets from slightly related distributions (40% edge agreement)
  set.seed(123)
  res.14 = perform_jointGHS_simulation(K.3,n.vals.3, p, N, frac.disagreement = fracs.disagreement[4], include.GHS=FALSE,include.SSJGL = F, include.JGL = F,nCores = nCores)
  
  # Case 15: datasets from unrelated distributions
  set.seed(123)
  res.15 = perform_jointGHS_simulation(K.3,n.vals.3, p, N, frac.disagreement = fracs.disagreement[5], include.GHS=FALSE,include.SSJGL = F, include.JGL = F,nCores = nCores)
  
  # Save results
  res.K10.list = list(res.11, res.12, res.13, res.14, res.15)  
  save(res.K10.list, file="data/jointGHS_simulations_K10.Rdata")
  
  # Print results
  #print_results_jointGHS(res.K10.list,fracs.disagreement, show.interval=F, show.sd=T,include.GHS=FALSE, include.SSJGL = F, include.JGL = F,collapse=T)
  
}



