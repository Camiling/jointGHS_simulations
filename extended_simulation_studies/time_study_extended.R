source('simulation_functions/perform_fastGHS_simulation.R')
source('simulation_functions/help_functions.R')
source('simulation_functions/perform_time_simulation.R')

# Write print function

nCores = 2 # If using HPC

# Perform simulation study

perform_time_sim_1 = FALSE
perform_time_sim_2 = FALSE
perform_time_sim_3 = FALSE
perform_time_sim_1_joint = FALSE
perform_time_sim_2_joint = FALSE
perform_time_sim_3_joint = FALSE
perform_time_sim_4_joint = TRUE
perform_time_sim_5_joint = FALSE
perform_time_sim_1_joint_smalln = FALSE
perform_time_sim_2_joint_smalln = FALSE
perform_time_sim_3_joint_smalln = FALSE


if(perform_time_sim_1){
  p=c(800)
  n=500
  time.res.large = perform_time_simulation(p,n,nCores, include.GHS=FALSE)
  save(time.res.large, file="extended_simulation_studies/data/time_simulations_extended_1.Rdata")
}

if(perform_time_sim_2){
  p=c(900)
  n=500
  time.res.large = perform_time_simulation(p,n,nCores, include.GHS=FALSE)
  save(time.res.large, file="extended_simulation_studies/data/time_simulations_extended_2.Rdata")
}

if(perform_time_sim_3){
  p=c(1000)
  n=500
  time.res.large = perform_time_simulation(p,n,nCores, include.GHS=FALSE)
  save(time.res.large, file="extended_simulation_studies/data/time_simulations_extended_3.Rdata")
}



# JOINT SIMULATION

if(perform_time_sim_1_joint){
  p=c(600,700)
  K=2
  n.vals=c(500, 500)
  time.res.large = perform_time_simulation_joint_onlyjointGHS(p,K=2,n.vals,nCores)
  save(time.res.large, file="extended_simulation_studies/data/time_simulations_extended_1_joint.Rdata")
}

if(perform_time_sim_2_joint){
  p=c(800,900)
  K=2
  n.vals=c(500,500)
  time.res.large = perform_time_simulation_joint_onlyjointGHS(p,K=2,n.vals,nCores)
  save(time.res.large, file="extended_simulation_studies/data/time_simulations_extended_2_joint.Rdata")
}

if(perform_time_sim_3_joint){
  p=c(1000)
  K=2
  n.vals=c(500, 500)
  time.res.large = perform_time_simulation_joint_onlyjointGHS(p,K=2,n.vals,nCores)
  save(time.res.large, file="extended_simulation_studies/data/time_simulations_extended_3_joint.Rdata")
}

if(perform_time_sim_4_joint){
  p=c(1100,1200)
  K=2
  n.vals=c(500, 500)
  time.res.large = perform_time_simulation_joint_onlyjointGHS(p,K=2,n.vals,nCores)
  save(time.res.large, file="extended_simulation_studies/data/time_simulations_extended_4_joint.Rdata")
}

if(perform_time_sim_5_joint){
  p=c(1500,2000)
  K=2
  n.vals=c(500, 500)
  time.res.large = perform_time_simulation_joint_onlyjointGHS(p,K=2,n.vals,nCores)
  save(time.res.large, file="extended_simulation_studies/data/time_simulations_extended_5_joint.Rdata")
}


# Smaller n

if(perform_time_sim_1_joint_smalln){
  p=c(800,900)
  K=2
  n.vals=c(200, 200)
  time.res.large = perform_time_simulation_joint_onlyjointGHS(p,K=2,n.vals,nCores)
  save(time.res.large, file="extended_simulation_studies/data/time_simulations_extended_1_joint_smalln.Rdata")
}

if(perform_time_sim_2_joint_smalln){
  p=c(1000,1100)
  K=2
  n.vals=c(200,200)
  time.res.large = perform_time_simulation_joint_onlyjointGHS(p,K=2,n.vals,nCores)
  save(time.res.large, file="extended_simulation_studies/data/time_simulations_extended_2_joint_smalln.Rdata")
}

if(perform_time_sim_3_joint_smalln){
  p=c(1200,1300)
  K=2
  n.vals=c(200, 200)
  time.res.large = perform_time_simulation_joint_onlyjointGHS(p,K=2,n.vals,nCores)
  save(time.res.large, file="extended_simulation_studies/data/time_simulations_extended_3_joint_smalln.Rdata")
}