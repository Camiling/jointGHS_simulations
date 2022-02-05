library(jointGHS) # Must be installed from Camiling/jointGHS on github
library(fastGHS) # Must be installed from Camiling/fastGHS on github
library(tailoredGlasso) # Must be installed from Camiling/tailoredGlasso on github
library(huge)
library(igraph)
library(JGL)
library(ggplot2)
library(gridExtra)
library(doParallel)
library(foreach)
library(parallel)
library(Rcpp)
source("SSJGL/R/SSJGL.R")
source("SSJGL/R/JGL.R")
source("SSJGL/R/admm.iters.R")
source("SSJGL/R/eval.R")
source("SSJGL/R/gete.R")

#' Perform jointGHS simulations
#' 
#' This function performs simulations for jointGHS, averaging over the results. 
#' 
#' @param K number of data sets/networks
#' @param n.vals the number of observations in each data set. A vector
#' @param p the number of nodes
#' @param N the number of simulations to perform
#' @param seeds seeds to use in each of the \eqn{N} simulations. A vector of length \eqn{N}. 
#' @param nCores how many cores should be used
#' @param frac.disagreement the fraction of edges that the networks should disagree on
#' @param method how should the similarity between the prec matrices be? Symmetric by default, meaning all are equally different. If not, one will stand out as completely unrelated.
#' @param include.jointGHS should the jointGHS be performed?
#' @param include.GHS should single-network GHS be performed?
#' @param include.SSJGL should we perform SSJGL?
#' @param include.JGL should the joint graphical lasso be included?
#' @param penalize.diagonal should the diagonal be penalized in the graphical lasso-based methods?
#' @param verbose logical indicator of printing information at each iteration
#' @param scale should the data be scaled?
#' @param save.res.jointGHS should Lambda, Nu and Theta be saved for all N simulations? Requires a lot of memory. Default FALSE.
#' @return simulation results, including sparsity, precision, recall and specificity
perform_jointGHS_simulation = function(K, n.vals, p, N=100, seeds=sample(1:1000,N), nCores = 3, frac.disagreement = 0, method='symmetric', 
                                       include.jointGHS=TRUE, include.GHS=TRUE, include.SSJGL=TRUE, include.JGL = TRUE, penalize.diagonal=FALSE ,verbose=TRUE, scale=TRUE, 
                                       save.res.jointGHS = FALSE){
 
  res=list()
  # JointGHS results
  res$opt.sparsities = matrix(0,N,K)
  res$precisions =  matrix(0,N,K)
  res$specificities =  matrix(0,N,K)
  res$recalls =  matrix(0,N,K)
  res$matrix.distances =  matrix(0,N,K)
  if(save.res.jointGHS){
    res$theta = list()
    res$E_NuInv = list()
    res$Lambda_sq = list()
  }
  # SSJGL results 
  res$opt.sparsities.ssjgl = matrix(0,N,K)
  res$precisions.ssjgl =  matrix(0,N,K)
  res$specificities.ssjgl =  matrix(0,N,K)
  res$recalls.ssjgl =  matrix(0,N,K)
  res$matrix.distances.ssjgl =  matrix(0,N,K)
  
  # JGL results (tuned by AIC)
  res$opt.sparsities.jgl = matrix(0,N,K)
  res$precisions.jgl =  matrix(0,N,K)
  res$specificities.jgl =  matrix(0,N,K)
  res$recalls.jgl =  matrix(0,N,K)
  res$matrix.distances.jgl =  matrix(0,N,K)
  
  # GHS results
  res$opt.sparsities.ghs = matrix(0,N,K)
  res$precisions.ghs =  matrix(0,N,K)
  res$specificities.ghs =  matrix(0,N,K)
  res$recalls.ghs =  matrix(0,N,K)
  res$matrix.distances.ghs =  matrix(0,N,K)
  
  # Start by generating the precision matrices 
  cov.matrices = list()
  prec.matrices = list()
  
  # Start by generating the first prec matrix
  huge.init = huge.generator(n.vals[1],p,graph='scale-free',verbose = F,v=1,u=0.01)
  theta.init = huge.init$omega
  theta.init[which(abs(theta.init)<1e-5,arr.ind=T)] = 0
  spars.init = huge.init$sparsity
  cov.matrices[[1]] = huge.init$sigma
  prec.matrices[[1]] = theta.init
  
  # Added this 
  if(scale) cov.matrices[[1]] = cov2cor(cov.matrices[[1]])
  if(scale) prec.matrices[[1]] = cov2cor(prec.matrices[[1]])
  
  # Avoid rounding errors leading to matrices not being symmetric
  if(!matrixcalc::is.symmetric.matrix(cov.matrices[[1]])){
    cov.matrices[[1]] = round(cov.matrices[[1]],8)
  }
  if(method=='symmetric'){
    for(k in 2:K){
      valid=F
      while(!valid){ # Ensure valid precision matrices
        huge.tmp = mutate.graph(huge.init,frac.disagreement,scale)
        cov.matrices[[k]] = huge.tmp$cov.mat
        prec.matrices[[k]] = huge.tmp$prec.mat
        # Avoid rounding errors leading to matrices not being symmetric
        if(!matrixcalc::is.symmetric.matrix(cov.matrices[[k]])){
          cov.matrices[[k]] = round(cov.matrices[[k]],8)
        }
        valid = matrixcalc::is.positive.definite(cov.matrices[[k]])
      }
    }
  }
  else{ # First K-1 graphs are similar
    for(k in 2:(K-1)){
      huge.tmp = mutate.graph(huge.init,frac.disagreement,scale)
      cov.matrices[[k]] = huge.tmp$cov.mat
      prec.matrices[[k]] = huge.tmp$prec.mat
      # Avoid rounding errors leading to matrices not being symmetric
      if(!matrixcalc::is.symmetric.matrix(cov.matrices[[k]])){
        cov.matrices[[k]] = round(cov.matrices[[k]],8)
      }
    }
    # Last graph is completely different
    huge.tmp = mutate.graph(huge.init,fraction = 1, scale)
    cov.matrices[[K]] = huge.tmp$cov.mat
    prec.matrices[[K]] = huge.tmp$prec.mat
    # Avoid rounding errors leading to matrices not being symmetric
    if(!matrixcalc::is.symmetric.matrix(cov.matrices[[K]])){
      cov.matrices[[K]] = round(cov.matrices[[K]],8)
    }
  }
  registerDoParallel(nCores)
  res.list = foreach (i=1:N) %dopar% {
    jointGHS_simulation_one_iteration(n.vals=n.vals,cov.matrices=cov.matrices,prec.matrices=prec.matrices,scale=scale,
                                     include.jointGHS=include.jointGHS, include.GHS=include.GHS, include.SSJGL=include.SSJGL, include.JGL=include.JGL, 
                                     penalize.diagonal=penalize.diagonal,seed=seeds[i], 
                                     save.res.jointGHS = save.res.jointGHS);
  }
  registerDoSEQ()
  
  # Save results from each replicate
  for(i in 1:N){
    est.tmp = res.list[[i]]
    
    # Results from jointGHS
    if(include.jointGHS){
      res$opt.sparsities[i,] = est.tmp$opt.sparsities 
      res$matrix.distances[i,] = est.tmp$matrix.distances
      res$precisions[i,] = est.tmp$precisions
      res$recalls[i,] = est.tmp$recalls
      res$specificities[i,] =  est.tmp$specificities
    }
    if(save.res.jointGHS){
      res$E_NuInv[[i]] = est.tmp$E_NuInv # A K by K matrix
      res$theta[[i]] = est.tmp$theta # A list of length K
      res$Lambda_sq[[i]] = est.tmp$Lambda_sq # A list of length K
    }
    # Results from jgl
    if(include.JGL){
      res$opt.sparsities.jgl[i,] = est.tmp$opt.sparsities.jgl
      res$matrix.distances.jgl[i,] = est.tmp$matrix.distances.jgl
      res$precisions.jgl[i,] = est.tmp$precisions.jgl
      res$recalls.jgl[i,] = est.tmp$recalls.jgl
      res$specificities.jgl[i,] =  est.tmp$specificities.jgl 
    }
    # Results from ssjgl
    if(include.SSJGL){
      res$opt.sparsities.ssjgl[i,] = est.tmp$opt.sparsities.ssjgl
      res$matrix.distances.ssjgl[i,] = est.tmp$matrix.distances.ssjgl
      res$precisions.ssjgl[i,] = est.tmp$precisions.ssjgl
      res$recalls.ssjgl[i,] = est.tmp$recalls.ssjgl
      res$specificities.ssjgl[i,] =  est.tmp$specificities.ssjgl 
    }
    # Results from ghs
    if(include.GHS){
      res$opt.sparsities.ghs[i,] = est.tmp$opt.sparsities.ghs
      res$matrix.distances.ghs[i,] = est.tmp$matrix.distances.ghs
      res$precisions.ghs[i,] = est.tmp$precisions.ghs
      res$recalls.ghs[i,] = est.tmp$recalls.ghs
      res$specificities.ghs[i,] =  est.tmp$specificities.ghs
    }
  }
  # Mean results from jointGHS
  if(include.jointGHS){
      res$mean.opt.sparsities= colMeans(res$opt.sparsities)
      res$mean.precisions= colMeans(res$precisions)
      res$mean.recalls = colMeans(res$recalls)
      res$mean.specificities =  colMeans(res$specificities)
      res$mean.matrix.distances = colMeans(res$matrix.distances)
  }
  # Mean results from JGL
  if(include.JGL){
    res$mean.opt.sparsities.jgl = colMeans(res$opt.sparsities.jgl)
    res$mean.precisions.jgl = colMeans(res$precisions.jgl)
    res$mean.recalls.jgl = colMeans(res$recalls.jgl)
    res$mean.specificities.jgl =  colMeans(res$specificities.jgl)
    res$mean.matrix.distances.jgl = colMeans(res$matrix.distances.jgl)
  }
  # Mean results from SSJGL
  if(include.SSJGL){
    res$mean.opt.sparsities.ssjgl = colMeans(res$opt.sparsities.ssjgl)
    res$mean.precisions.ssjgl = colMeans(res$precisions.ssjgl)
    res$mean.recalls.ssjgl = colMeans(res$recalls.ssjgl)
    res$mean.specificities.ssjgl =  colMeans(res$specificities.ssjgl)
    res$mean.matrix.distances.ssjgl = colMeans(res$matrix.distances.ssjgl)
  }
  # Mean results from GHS
  if(include.GHS){
    res$mean.opt.sparsities.ghs = colMeans(res$opt.sparsities.ghs)
    res$mean.precisions.ghs = colMeans(res$precisions.ghs)
    res$mean.recalls.ghs = colMeans(res$recalls.ghs)
    res$mean.specificities.ghs =  colMeans(res$specificities.ghs)
    res$mean.matrix.distances.ghs = colMeans(res$matrix.distances.ghs)
  }
  
  res$true.prec.matrices = prec.matrices
  res$true.sparsity = spars.init
  return(res)
  
}

# Function for performing one iteration -----------------------------------------

jointGHS_simulation_one_iteration = function(n.vals,cov.matrices,prec.matrices,scale,include.jointGHS, include.GHS, include.SSJGL,include.JGL, 
                                            penalize.diagonal,seed, save.res.jointGHS) {
  y = list()
  K=length(n.vals)
  p=ncol(prec.matrices[[1]])
  ghs.res = list()
  set.seed(seed)
  # Generate data. 
  for(k in 1:K){
    y[[k]] = mvtnorm::rmvnorm(n.vals[k], mean=rep(0,p), cov.matrices[[k]])
    if (scale) y[[k]] = scale(y[[k]])
    # Use single-network method
    if(include.GHS){
      ghs.res[[k]] = fastGHS::fastGHS(X=y[[k]], AIC_selection = T, AIC_eps = 1, epsilon = 1e-3, verbose = F)$theta
    }
  }
  # Perform joint methods
  if(include.jointGHS){
    res.full.tmp = jointGHS(X=y, AIC_selection = T, AIC_eps = 0.1, epsilon=1e-3, verbose = F)
    res.tmp = res.full.tmp$theta
    #if(include.GHS){
    #  ghs.res = res.full.tmp$theta_single
    #}
  }
  if(include.SSJGL){
    lam1 = 1
    lam2 = 1
    v1 = 1
    lam.eff = lam1 + c(1:10) * 5
    v0s = lam1/lam.eff
    ssjgl.tmp = SSJGL(Y=y,penalty='fused',lambda0=1, lambda1=lam1,lambda2=lam2, v1 = v1, v0s = v0s, tol.em=1e-4, a=1, b=p, doubly=TRUE, normalize=TRUE)
    ssjgl.tmp = ssjgl.tmp$thetalist[[10]]
  }
  if(include.JGL){
    jgl.tmp = JGL_select_AIC(Y=y,penalty='fused',nlambda1=10,lambda1.min=0.01,lambda1.max=1,nlambda2=10,lambda2.min=0,lambda2.max=0.1,lambda2.init=0.01,
                             penalize.diagonal=penalize.diagonal)
  }
  est=list()
  
  # If save.res.jointGHS==T, save Lambda, expectation of 1/Nu and Theta
  if(save.res.jointGHS){
    est$E_NuInv = res.full.tmp$E_NuInv # A K by K matrix
    est$theta = res.full.tmp$theta # A list of length K
    est$Lambda_sq = res.full.tmp$Lambda_sq # A list of length K
  }
  
  # Results from jointGHS
  if(include.jointGHS){
    res.tmp = lapply(res.tmp, cov2cor)
    #est$opt.sparsities  = unlist(lapply(res.tmp, FUN= function(m) tailoredGlasso::sparsity(abs(res.tmp[[k]])>1e-5)))
    est$opt.sparsities  = unlist(lapply(res.tmp, FUN= function(m) tailoredGlasso::sparsity(abs(m)>1e-5)))
    est$matrix.distances = sapply(1:K,FUN=function(k) matrix.distance(prec.matrices[[k]], res.tmp[[k]]))
    est$precisions =  sapply(1:K,FUN=function(k) precision(prec.matrices[[k]]!=0, abs(res.tmp[[k]])>1e-5))
    est$recalls =  sapply(1:K,FUN=function(k) recall(prec.matrices[[k]]!=0, abs(res.tmp[[k]])>1e-5 ))
    est$specificities =  sapply(1:K,FUN=function(k) specificity(prec.matrices[[k]]!=0, abs(res.tmp[[k]])>1e-5 ))
  }
  # Results from SSJGL
  if(include.SSJGL){
    ssjgl.tmp = lapply(ssjgl.tmp, cov2cor)
    est$opt.sparsities.ssjgl  = unlist(lapply(ssjgl.tmp, FUN= function(m) tailoredGlasso::sparsity(abs(m)>1e-5)))
    est$matrix.distances.ssjgl = sapply(1:K,FUN=function(k) matrix.distance(prec.matrices[[k]], ssjgl.tmp[[k]]))
    est$precisions.ssjgl = sapply(1:K,FUN=function(k) precision(prec.matrices[[k]]!=0, ssjgl.tmp[[k]]!=0))
    est$recalls.ssjgl = sapply(1:K,FUN=function(k) recall(prec.matrices[[k]]!=0, ssjgl.tmp[[k]]!=0))
    est$specificities.ssjgl = sapply(1:K,FUN=function(k) specificity(prec.matrices[[k]]!=0, ssjgl.tmp[[k]]!=0))
  }
  # Results from JGL
  if(include.JGL){
    est$opt.sparsities.jgl = jgl.tmp$opt.sparsities
    est$matrix.distances.jgl = sapply(1:K,FUN=function(k) matrix.distance(prec.matrices[[k]], jgl.tmp$opt.fit[[k]]))
    est$precisions.jgl = sapply(1:K,FUN=function(k) precision(prec.matrices[[k]]!=0, jgl.tmp$opt.fit[[k]]!=0))
    est$recalls.jgl = sapply(1:K,FUN=function(k) recall(prec.matrices[[k]]!=0, jgl.tmp$opt.fit[[k]]!=0))
    est$specificities.jgl = sapply(1:K,FUN=function(k) specificity(prec.matrices[[k]]!=0, jgl.tmp$opt.fit[[k]]!=0))
  }
  # Results from GHS
  if(include.GHS){
    est$opt.sparsities.ghs = unlist(lapply(ghs.res, FUN= function(m) tailoredGlasso::sparsity(abs(m)>1e-5)))
    est$matrix.distances.ghs = sapply(1:K,FUN=function(k) matrix.distance(prec.matrices[[k]], ghs.res[[k]]))
    est$precisions.ghs = sapply(1:K,FUN=function(k) precision(prec.matrices[[k]]!=0, abs(ghs.res[[k]])>1e-5))
    est$recalls.ghs = sapply(1:K,FUN=function(k) recall(prec.matrices[[k]]!=0, abs(ghs.res[[k]])>1e-5))
    est$specificities.ghs =  sapply(1:K,FUN=function(k) specificity(prec.matrices[[k]]!=0, abs(ghs.res[[k]])>1e-5))
  }
  return(est)
}




