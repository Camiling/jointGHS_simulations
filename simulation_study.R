source('simulation_functions/perform_jointGHS_simulation.R')
source('simulation_functions/help_functions.R')

nCores = 56 # If ran on HPC
registerDoParallel(nCores) # Moved to outside of function, after Colin's suggestion

# Perform simulation study

# 100 simulations per case
N = 10
p=100

# K = 2 data sets, of various similarity -------------
K=2
n.vals = c(150,200)
#n.vals = c(100,150)
fracs.disagreement = c(0,0.2,0.4,1)

# Case 1: datasets from same distribution
set.seed(1234)
res.1 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = fracs.disagreement[1], tau_sq = 10, tau_sq_ghs = 0.01, ebic.gamma=0,method='symmetric')

# Case 2: datasets from similar distributions (80% edge agreement)
set.seed(1234)
res.2 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = fracs.disagreement[2], tau_sq = 10, tau_sq_ghs = 0.01, ebic.gamma=0,method='symmetric')

# Case 3: datasets from slightly related distributions (60% edge agreement)
set.seed(1234)
res.3 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = fracs.disagreement[3], tau_sq = 10, tau_sq_ghs = 0.01, ebic.gamma=0,method='symmetric')

# Case 4: datasets from unrelated distributions
set.seed(1234)
res.4 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = fracs.disagreement[4], tau_sq = 10, tau_sq_ghs = 0.01, ebic.gamma=0,method='symmetric')

# Print results
print_results_jointGHS_show_SD(list(res.1, res.2, res.3, res.4),fracs.disagreement)

save(list(res.1, res.2, res.3, res.4), 'jointGHS_simulations_K2.Rdata')

# K = 4 data sets, of various similarity -------------
K.2 = 4
n.vals.2 = c(150, 200, 150, 100)

# Case 5: datasets from same distribution
set.seed(123)
res.5 = perform_jointGHS_simulation(K.2,n.vals.2, p, N, frac.disagreement = fracs.disagreement[1], tau_sq = 100, tau_sq_ghs = 0.01, ebic.gamma=0,method='symmetric')

# Case 6: datasets from similar distributions (80% edge agreement)
set.seed(123)
res.6 = perform_jointGHS_simulation(K.2,n.vals.2, p, N, frac.disagreement = fracs.disagreement[2], tau_sq = 100, tau_sq_ghs = 0.01, ebic.gamma=0,method='symmetric')

# Case 7: datasets from slightly related distributions (60% edge agreement)
set.seed(123)
res.7 = perform_jointGHS_simulation(K.2,n.vals.2, p, N, frac.disagreement = fracs.disagreement[3], tau_sq = 100, tau_sq_ghs = 0.01, ebic.gamma=0,method='symmetric')

# Case 8: datasets from unrelated distributions
set.seed(123)
res.8 = perform_jointGHS_simulation(K.2,n.vals.2, p, N, frac.disagreement = fracs.disagreement[4], tau_sq = 100, tau_sq_ghs = 0.01, ebic.gamma=0,method='symmetric')

# Print results
print_results_jointGHS_show_SD(list(res.5, res.6, res.7, res.8),fracs.disagreement)

save(list(res.5, res.6, res.7, res.8), 'jointGHS_simulations_K4.Rdata')

# K = 10 data sets, of various similarity -------------
K.3=10
n.vals.3 = c(150,200, 150, 100, 120, 140, 170, 180, 200, 190)

# Case 9: datasets from same distribution
set.seed(123)
res.9 = perform_jointGHS_simulation(K.3,n.vals.3, p, N, frac.disagreement = fracs.disagreement[1], tau_sq = 100, tau_sq_ghs = 0.01, ebic.gamma=0,method='symmetric')

# Case 10: datasets from similar distributions (80% edge agreement)
set.seed(123)
res.10 = perform_jointGHS_simulation(K.3,n.vals.3, p, N, frac.disagreement = fracs.disagreement[2], tau_sq = 100, tau_sq_ghs = 0.01, ebic.gamma=0,method='symmetric')

# Case 11: datasets from slightly related distributions (60% edge agreement)
set.seed(123)
res.11 = perform_jointGHS_simulation(K.3,n.vals.3, p, N, frac.disagreement = fracs.disagreement[3], tau_sq = 100, tau_sq_ghs = 0.01, ebic.gamma=0,method='symmetric')

# Case 12: datasets from unrelated distributions
set.seed(123)
res.12 = perform_jointGHS_simulation(K.3,n.vals.3, p, N, frac.disagreement = fracs.disagreement[4], tau_sq = 100, tau_sq_ghs = 0.01, ebic.gamma=0,method='symmetric')

# Print results
print_results_jointGHS_show_SD(list(res.9, res.10, res.11, res.12),fracs.disagreement)

save(list(res.9, res.10, res.11, res.12), 'jointGHS_simulations_K10.Rdata')

registerDoSEQ()

# Save all simulation results
res.list = list(res.1, res.2, res.3, res.4, res.5, res.6, res.7, res.8, res.9, res.10, res.11, res.12)
save(res.list, 'jointGHS_simulations.Rdata')








