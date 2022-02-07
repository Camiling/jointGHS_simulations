source('simulation_functions/perform_fastGHS_simulation.R')
source('simulation_functions/help_functions.R')
source('simulation_functions/perform_time_simulation.R')

# Plot results from time simulations

load("data/time_simulations_small.Rdata")
load("data/time_simulations_large.Rdata")

# GHS vs fastGHS

df.time = cbind(time.res[[1]], time.res[[2]])