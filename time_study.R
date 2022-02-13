source('simulation_functions/perform_fastGHS_simulation.R')
source('simulation_functions/help_functions.R')
source('simulation_functions/perform_time_simulation.R')

# Write print function

nCores = 8 # If using HPC

# Perform simulation study

perform_time_sim_small = FALSE
perform_time_sim_large= TRUE


if(perform_time_sim_small){
  p=c(10,20,30,40,50,60,70,80,90)
  n=100
  time.res = perform_time_simulation(p,n,nCores, AIC_selection=FALSE)
  save(time.res, file="data/time_simulations_small.Rdata")
}

if(perform_time_sim_large){
  #p=c(20,50,100,200,300,400,500,600,700,800,900,1000)
  p=c(20,50,100,200,300,400,500,600)
  n=500
  time.res.large = perform_time_simulation(p,n,nCores, include.GHS=FALSE)
  save(time.res.large, file="data/time_simulations_large.Rdata")
}


