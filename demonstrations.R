library(jointGHS) # Must be installed from Camiling/jointGHS on github
library(fastGHS) # Must be installed from Camiling/fastGHS
library(tailoredGlasso) # Must be installed from Camiling/tailoredGlasso
library(huge)
library(glasso)
library(igraph)
library(JGL)
library(ggplot2)
library(gridExtra)
library(doParallel)
library(foreach)
library(parallel)
library(Rcpp)
source('simulation_functions/perform_jointGHS_simulation.R')
source('simulation_functions/help_functions.R')

# Demonstrate that JGL is not computationally feasible for K=10 -----------------------------

set.seed(123456)
K=10
p=100
n.vals = c(150,200, 130, 160, 180, 170, 200, 150, 190, 150)
cov.matrices = list()
prec.matrices = list()
# Start by generating the first prec matrix
huge.init = huge.generator(n.vals[1],p,graph='scale-free',verbose = F,v=0.5,u=0.05)
theta.init = huge.init$omega
theta.init[which(abs(theta.init)<1e-5,arr.ind=T)] = 0
spars.init = huge.init$sparsity
cov.matrices[[1]] = huge.init$sigma
prec.matrices[[1]] = theta.init
# Added this 
cov.matrices[[1]] = cov2cor(cov.matrices[[1]])
# Avoid rounding errors leading to matrices not being symmetric
if(!matrixcalc::is.symmetric.matrix(cov.matrices[[1]])){
  cov.matrices[[1]] = round(cov.matrices[[1]],8)
}
for(k in 2:K){
  huge.tmp = mutate.graph(huge.init,0.1, scale=F)
  cov.matrices[[k]] = huge.tmp$cov.mat
  prec.matrices[[k]] = huge.tmp$prec.mat
  # Avoid rounding errors leading to matrices not being symmetric
  if(!matrixcalc::is.symmetric.matrix(cov.matrices[[k]])){
    cov.matrices[[k]] = round(cov.matrices[[k]],8)
  }
}
y = list()
# Generate data. 
scale=T
for(k in 1:K){
  y[[k]] = mvtnorm::rmvnorm(n.vals[k], mean=rep(0,p), cov.matrices[[k]])
  if (scale) y[[k]] = scale(y[[k]])
}


jgl.tmp = JGL_select_AIC(Y=y,penalty='fused',nlambda1=10,lambda1.min=0.01,lambda1.max=1,nlambda2=10,lambda2.min=0,lambda2.max=0.1,lambda2.init=0.01,
                         penalize.diagonal=F)

# Still running at >20 hrs => N=100 simulation case would take more than 80 days if it even is feasible at all. And this is just for a 10x10 grid of lambda1 and lambda2



# Demonstrate that SSL does not work for K=10 ---------------------------------

source("SSJGL/R/SSJGL.R")
source("SSJGL/R/JGL.R")
source("SSJGL/R/admm.iters.R")
source("SSJGL/R/eval.R")
source("SSJGL/R/gete.R")

penalty <- "fused"
lambda1 <- 1 
lambda2 <- 1 
v1 <- 1
lambda.eff2 <- lambda1 + 1 + c(0:20)*10
v0s <- lambda1/lambda.eff2
fit.ssjgl = SSJGL(Y=y,penalty=penalty,lambda0=1, lambda1=lambda1,lambda2=lambda2, v1 = v1, v0s = v0s, tol.em=0.0001, a=1, b=1, doubly=FALSE, c = 0.01)

# Also infeasbile, which makes sense given that it is an adaption of the JGL

# What about K=4?

y.2 = list(y[[1]], y[[2]], y[[3]], y[[4]])
fit.ssjgl = SSJGL(Y=y.2,penalty=penalty,lambda0=1, lambda1=lambda1,lambda2=lambda2, v1 = v1, v0s = v0s, tol.em=0.0001, a=1, b=1, doubly=FALSE, c = 0.01)

# This is feasible.


