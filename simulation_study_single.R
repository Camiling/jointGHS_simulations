source('simulation_functions/perform_fastGHS_simulation.R')
source('simulation_functions/help_functions.R')

# Write this study

# Write print function

nCores = 5 # If using HPC

# Perform simulation study


# 100 simulations per case
N = 10

# Case 1: p=50,n=100
set.seed(12345)
p = 50
n = 100
res.1 = perform_fastGHS_simulation(n, p, N,  nCores = nCores)
  
# Case 2: p=200,n=150 (Here GHS with Gibbs fails, so we do not run)
set.seed(12345)
p = 200
n = 150
res.2 = perform_fastGHS_simulation(n, p, N, include.GHS = F, nCores = nCores)
  
# Case 3: p=400, n=300
set.seed(1234)
p = 400
n = 300
res.3 = perform_fastGHS_simulation(n, p, N, nCores = nCores, include.GHS = F)
  
# Case 4: p=100, n=100
set.seed(1234)
p = 100
n = 100
res.4 = perform_fastGHS_simulation(n, p, N, nCores = nCores, include.GHS = F)

# Case 5: p=50, n=200
set.seed(1234)
p = 50
n = w00
res.4 = perform_fastGHS_simulation(n, p, N, nCores = nCores, include.GHS = F)
  
# Save results
res.fast = list(res.1, res.2, res.3, res.4, res.5)
save(res.fast, file="data/fastGHS_simulations.Rdata")
  
# Print results
#print_results_fastGHS(res.fast, show.interval=F, show.sd=F)