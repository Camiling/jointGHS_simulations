source('simulation_functions/perform_jointGHS_simulation.R')
source('simulation_functions/help_functions.R')

nCores = 56 # In using HPC

# Perform simulation study

# Choose which scenarios to simulate
run.K2 = F
run.K4 = F
run.K10 = T
run.K4.oneout = F

# 100 simulations per case
N = 100
p = 100

# Set parameters
stars.thresh = 0.01 # Used in the ordinary graphical lasso
ebic.gamma = 0.2 # Used in JoStARS
var.thresh.jostars = 0.05 # Used in JoStARS

# K = 2 data sets, of various similarity -------------
K=2
n.vals = c(150,200)
fracs.disagreement = c(0,0.2,0.4,1)

if(run.K2){
  # Case 1: datasets from same distribution
  set.seed(1234)
  res.1 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = fracs.disagreement[1], tau_sq = 10, tau_sq_ghs = 0.01,ebic.gamma=ebic.gamma,stars.thresh=stars.thresh, var.thresh.jostars = var.thresh.jostars, method='symmetric', nCores = nCores)
  
  # Case 2: datasets from similar distributions (80% edge agreement)
  set.seed(1234)
  res.2 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = fracs.disagreement[2], tau_sq = 10, tau_sq_ghs = 0.01,  ebic.gamma=ebic.gamma,stars.thresh=stars.thresh, var.thresh.jostars = var.thresh.jostars,method='symmetric', nCores = nCores)
  
  # Case 3: datasets from slightly related distributions (60% edge agreement)
  set.seed(1234)
  res.3 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = fracs.disagreement[3], tau_sq = 10, tau_sq_ghs = 0.01,ebic.gamma=ebic.gamma,stars.thresh=stars.thresh, var.thresh.jostars = var.thresh.jostars,method='symmetric', nCores = nCores)
  
  # Case 4: datasets from unrelated distributions
  set.seed(1234)
  res.4 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = fracs.disagreement[4], tau_sq = 10, tau_sq_ghs = 0.01, ebic.gamma=ebic.gamma,stars.thresh=stars.thresh, var.thresh.jostars = var.thresh.jostars,method='symmetric', nCores = nCores)
  
  # Save results
  res.K2 = list(res.1, res.2, res.3, res.4)
  save(res.K2, file="data/jointGHS_simulations_K2.Rdata")
  
  # Print results
  print_results_jointGHS(res.K2, fracs.disagreement, include.jostars = TRUE, show.interval=F, show.sd=T)
}




# K = 4 data sets, of various similarity -------------
K.2 = 4
n.vals.2 = c(150, 200, 150, 100)

# Went from 1234 seed to 12345

if(run.K4){
  # Case 5: datasets from same distribution
  set.seed(12345)
  res.5 = perform_jointGHS_simulation(K.2,n.vals.2, p, N, frac.disagreement = fracs.disagreement[1], tau_sq = 10, tau_sq_ghs = 0.008, ebic.gamma=ebic.gamma,stars.thresh=stars.thresh, var.thresh.jostars = var.thresh.jostars, method='symmetric', nCores = nCores)
  
  # Case 6: datasets from similar distributions (80% edge agreement)
  set.seed(12345)
  res.6 = perform_jointGHS_simulation(K.2,n.vals.2, p, N, frac.disagreement = fracs.disagreement[2], tau_sq = 10, tau_sq_ghs = 0.008,  ebic.gamma=ebic.gamma,stars.thresh=stars.thresh, var.thresh.jostars = var.thresh.jostars, method='symmetric', nCores = nCores)
  
  # Case 7: datasets from slightly related distributions (60% edge agreement)
  set.seed(12345)
  res.7 = perform_jointGHS_simulation(K.2,n.vals.2, p, N, frac.disagreement = fracs.disagreement[3], tau_sq = 10, tau_sq_ghs = 0.008,ebic.gamma=ebic.gamma,stars.thresh=stars.thresh, var.thresh.jostars = var.thresh.jostars, method='symmetric', nCores = nCores)
  
  # Case 8: datasets from unrelated distributions
  set.seed(12345)
  res.8 = perform_jointGHS_simulation(K.2,n.vals.2, p, N, frac.disagreement = fracs.disagreement[4], tau_sq = 10, tau_sq_ghs = 0.008, ebic.gamma=ebic.gamma,stars.thresh=stars.thresh, var.thresh.jostars = var.thresh.jostars, method='symmetric', nCores = nCores)
  
  # Save results
  res.K4.list = list(res.5, res.6, res.7, res.8)
  save(res.K4.list, file="data/jointGHS_simulations_K4.Rdata")
  
  # Print results
  print_results_jointGHS(res.K4.list,fracs.disagreement, include.jostars = F, show.interval=F, show.sd=T)
}


# K = 10 data sets, of various similarity -------------
K.3=10
#n.vals.3 = c(150,200, 150, 100, 120, 140, 170, 180, 200, 190)
n.vals.3 = c(150,100, 130, 100, 120, 140, 150, 110, 120, 100)


if(run.K10){
  # Case 9: datasets from same distribution
  set.seed(123)
  res.9 = perform_jointGHS_simulation(K.3,n.vals.3, p, N, frac.disagreement = fracs.disagreement[1], tau_sq = 100, tau_sq_ghs = 0.01, stars.thresh=stars.thresh, include.jostars = FALSE, method='symmetric', nCores = nCores)
  
  # Case 10: datasets from similar distributions (80% edge agreement)
  set.seed(123)
  res.10 = perform_jointGHS_simulation(K.3,n.vals.3, p, N, frac.disagreement = fracs.disagreement[2], tau_sq = 100, tau_sq_ghs = 0.01, stars.thresh=stars.thresh, include.jostars = FALSE, method='symmetric', nCores = nCores)
  
  # Case 11: datasets from slightly related distributions (60% edge agreement)
  set.seed(123)
  res.11 = perform_jointGHS_simulation(K.3,n.vals.3, p, N, frac.disagreement = fracs.disagreement[3], tau_sq = 100, tau_sq_ghs = 0.01,stars.thresh=stars.thresh, include.jostars = FALSE, method='symmetric', nCores = nCores)
  
  # Case 12: datasets from unrelated distributions
  set.seed(123)
  res.12 = perform_jointGHS_simulation(K.3,n.vals.3, p, N, frac.disagreement = fracs.disagreement[4], tau_sq = 100, tau_sq_ghs = 0.01,stars.thresh=stars.thresh, include.jostars = FALSE, method='symmetric', nCores = nCores)
  
  # Save results
  res.K10.list = list(res.9, res.10, res.11, res.12)  
  save(res.K10.list, file="data/jointGHS_simulations_K10.Rdata")
  
  # Print results
  print_results_jointGHS(res.K10.list,fracs.disagreement, include.jostars=F, show.interval=F, show.sd=T, collapse=T)

}


# K = 4 data sets, of various similarity, with one standing out as unrelated -------------
K.4 = 4
n.vals.4 = c(150, 200, 150, 200)

if(run.K4.oneout){
  # Case 13: datasets from same distribution
  set.seed(12345)
  res.13 = perform_jointGHS_simulation(K.4,n.vals.4, p, N, frac.disagreement = fracs.disagreement[1], tau_sq = 10, tau_sq_ghs = 0.01, ebic.gamma=ebic.gamma,stars.thresh=stars.thresh, var.thresh.jostars = var.thresh.jostars, method='unsymmetric', nCores = nCores)
  
  # Case 14: datasets from similar distributions (80% edge agreement)
  set.seed(12345)
  res.14 = perform_jointGHS_simulation(K.4,n.vals.4, p, N, frac.disagreement = fracs.disagreement[2], tau_sq = 10, tau_sq_ghs = 0.01,  ebic.gamma=ebic.gamma,stars.thresh=stars.thresh, var.thresh.jostars = var.thresh.jostars, method='unsymmetric', nCores = nCores)
  
  # Case 15: datasets from slightly related distributions (60% edge agreement)
  set.seed(12345)
  res.15 = perform_jointGHS_simulation(K.4,n.vals.4, p, N, frac.disagreement = fracs.disagreement[3], tau_sq = 10, tau_sq_ghs = 0.01,ebic.gamma=ebic.gamma,stars.thresh=stars.thresh, var.thresh.jostars = var.thresh.jostars, method='unsymmetric', nCores = nCores)
  
  # Save results
  res.K4.oneout.list = list(res.13, res.14, res.15)
  save(res.K4.oneout.list, file="data/jointGHS_simulations_K4_oneout.Rdata")
  
  # Print results
  print_results_jointGHS(res.K4.oneout.list,fracs.disagreement, include.jostars = TRUE, show.interval=F, show.sd=T)
}





