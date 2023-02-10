rm(list=ls())
library(Rcpp)
library(JGL)
library(huge)
library(igraph)
library(ggplot2)
library(gridExtra)
library(doParallel)
library(foreach)
library(parallel)
library(Rcpp)
library(tailoredGlasso) 
source('simulation_functions/help_functions.R')
source('simulation_functions/perform_time_simulation.R')
source('GemBag/Gembag_implementation/BIC_GemBag.R')
Rcpp::sourceCpp('GemBag/Gembag_implementation/GemBag-algo.cpp')
source("SSJGL/R/SSJGL.R")
source("SSJGL/R/JGL.R")
source("SSJGL/R/admm.iters.R")
source("SSJGL/R/eval.R")
source("SSJGL/R/gete.R")


# Write print function

nCores = 6 # If using HPC

# Perform simulation study

perform_time_sim_joint= TRUE


if(perform_time_sim_joint){
  p=c(20,50,100,150,200,250)
  # First for K=2 networks
  K=2
  n.vals = c(100,150)
  set.seed(33)
  time.res.joint.K2= perform_time_simulation_joint(p,K,n.vals,nCores)
  save(time.res.joint.K2, file="extended_simulation_studies/data/time_simulations_joint_K2.Rdata")
  # First for K=3 networks
  K=3
  n.vals = c(100,100,150)
  set.seed(1234)
  time.res.joint.K3= perform_time_simulation_joint(p,K,n.vals,nCores)
  save(time.res.joint.K3, file="extended_simulation_studies/data/time_simulations_joint_K3.Rdata")
  # First for K=4 networks
  K=4
  n.vals = c(100,100,150,150)
  set.seed(123)
  time.res.joint.K4= perform_time_simulation_joint(p,K,n.vals,nCores)
  save(time.res.joint.K4, file="extended_simulation_studies/data/time_simulations_joint_K4.Rdata")
}


