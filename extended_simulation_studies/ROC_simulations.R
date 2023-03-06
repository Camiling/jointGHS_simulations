rm(list=ls())
source('extended_simulation_studies/perform_ROC_simulation.R')
source('simulation_functions/help_functions.R')

nCores = 25 # If using HPC
# Perform simulation study

# p=50 -----------------------------------------------------

# Choose which scenarios to simulate
run.jointGHS = F
run.JGL = F
run.SSJGL = F
run.GemBag = T
# Gembag: can just use the existing inclusion probs found in original sim.

# 100 simulations per case
N = 100
p = 50
K=2
n.vals = c(50,80)
fracs.disagreement = c(0,0.2,0.4,0.6,0.8,1)

if(run.jointGHS){
  set.seed(1234)
  res.1 = perform_jointGHS_simulation_ROC(K,n.vals, p, N, frac.disagreement = fracs.disagreement[1],method='jointGHS', nCores = nCores)
  set.seed(1234)
  res.2 = perform_jointGHS_simulation_ROC(K,n.vals, p, N, frac.disagreement = fracs.disagreement[2],method='jointGHS', nCores = nCores)
  set.seed(1234)
  res.3 = perform_jointGHS_simulation_ROC(K,n.vals, p, N, frac.disagreement = fracs.disagreement[3],method='jointGHS', nCores = nCores)
  set.seed(1234)
  res.4 = perform_jointGHS_simulation_ROC(K,n.vals, p, N, frac.disagreement = fracs.disagreement[4],method='jointGHS', nCores = nCores)
  set.seed(1234)
  res.5 = perform_jointGHS_simulation_ROC(K,n.vals, p, N, frac.disagreement = fracs.disagreement[5],method='jointGHS', nCores = nCores)
  set.seed(1234)
  res.6 = perform_jointGHS_simulation_ROC(K,n.vals, p, N, frac.disagreement = fracs.disagreement[6],method='jointGHS', nCores = nCores)
  # Save results
  res.ROC.jointGHS = list(res.1, res.2, res.3, res.4, res.5, res.6)
  save(res.ROC.jointGHS, file="extended_simulation_studies/data/jointGHS_simulations_jointGHS_ROC.Rdata")
}
  
if(run.JGL){
  set.seed(1234)
  res.1 = perform_jointGHS_simulation_ROC(K,n.vals, p, N, frac.disagreement = fracs.disagreement[1],method='JGL', nCores = nCores)
  set.seed(1234)
  res.2 = perform_jointGHS_simulation_ROC(K,n.vals, p, N, frac.disagreement = fracs.disagreement[2],method='JGL', nCores = nCores)
  set.seed(1234)
  res.3 = perform_jointGHS_simulation_ROC(K,n.vals, p, N, frac.disagreement = fracs.disagreement[3],method='JGL', nCores = nCores)
  set.seed(1234)
  res.4 = perform_jointGHS_simulation_ROC(K,n.vals, p, N, frac.disagreement = fracs.disagreement[4],method='JGL', nCores = nCores)
  set.seed(1234)
  res.5 = perform_jointGHS_simulation_ROC(K,n.vals, p, N, frac.disagreement = fracs.disagreement[5],method='JGL', nCores = nCores)
  set.seed(1234)
  res.6 = perform_jointGHS_simulation_ROC(K,n.vals, p, N, frac.disagreement = fracs.disagreement[6],method='JGL', nCores = nCores)
  # Save results
  res.ROC.JGL = list(res.1, res.2, res.3, res.4, res.5, res.6)
  save(res.ROC.JGL, file="extended_simulation_studies/data/jointGHS_simulations_JGL_ROC.Rdata")
}
  
if(run.SSJGL){
  set.seed(1234)
  res.1 = perform_jointGHS_simulation_ROC(K,n.vals, p, N, frac.disagreement = fracs.disagreement[1],method='SSJGL', nCores = nCores)
  set.seed(1234)
  res.2 = perform_jointGHS_simulation_ROC(K,n.vals, p, N, frac.disagreement = fracs.disagreement[2],method='SSJGL', nCores = nCores)
  set.seed(1234)
  res.3 = perform_jointGHS_simulation_ROC(K,n.vals, p, N, frac.disagreement = fracs.disagreement[3],method='SSJGL', nCores = nCores)
  set.seed(1234)
  res.4 = perform_jointGHS_simulation_ROC(K,n.vals, p, N, frac.disagreement = fracs.disagreement[4],method='SSJGL', nCores = nCores)
  set.seed(1234)
  res.5 = perform_jointGHS_simulation_ROC(K,n.vals, p, N, frac.disagreement = fracs.disagreement[5],method='SSJGL', nCores = nCores)
  set.seed(1234)
  res.6 = perform_jointGHS_simulation_ROC(K,n.vals, p, N, frac.disagreement = fracs.disagreement[6],method='SSJGL', nCores = nCores)
  # Save results
  res.ROC.SSJGL = list(res.1, res.2, res.3, res.4, res.5, res.6)
  save(res.ROC.SSJGL, file="extended_simulation_studies/data/jointGHS_simulations_SSJGL_ROC.Rdata")
}  

# Finally, the same for GemBag
if(run.GemBag){
  set.seed(1234)
  res.1 = perform_jointGHS_simulation_ROC(K,n.vals, p, N, frac.disagreement = fracs.disagreement[1],method='GemBag', nCores = nCores)
  set.seed(1234)
  res.2 = perform_jointGHS_simulation_ROC(K,n.vals, p, N, frac.disagreement = fracs.disagreement[2],method='GemBag', nCores = nCores)
  set.seed(1234)
  res.3 = perform_jointGHS_simulation_ROC(K,n.vals, p, N, frac.disagreement = fracs.disagreement[3],method='GemBag', nCores = nCores)
  set.seed(1234)
  res.4 = perform_jointGHS_simulation_ROC(K,n.vals, p, N, frac.disagreement = fracs.disagreement[4],method='GemBag', nCores = nCores)
  set.seed(1234)
  res.5 = perform_jointGHS_simulation_ROC(K,n.vals, p, N, frac.disagreement = fracs.disagreement[5],method='GemBag', nCores = nCores)
  set.seed(1234)
  res.6 = perform_jointGHS_simulation_ROC(K,n.vals, p, N, frac.disagreement = fracs.disagreement[6],method='GemBag', nCores = nCores)
  # Save results
  res.ROC.gembag = list(res.1, res.2, res.3, res.4, res.5, res.6)
  save(res.ROC.gembag, file="extended_simulation_studies/data/jointGHS_simulations_GemBag_ROC.Rdata")
} 

