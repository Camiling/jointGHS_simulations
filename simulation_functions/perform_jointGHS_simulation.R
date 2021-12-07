library(jointGHS) # Must be installed from Camiling/jointGHS on github
library(fastGHS) # Must be installed from Camiling/fastGHS
library(tailoredGlasso) # Must be installed from Camiling/tailoredGlasso
library(JoStARS) # Must be installed from Camiling/JoStARS
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
#' @param tau_sq the values of tau to use in the jointGHS. Either a vector of length \eqn{K} or a single value to use for all networks
#' @param tau_sq_ghs the value of tau to use in the GHS
#' @param ebic.gamma the additional penalty term in the extended BIC criterion for selecting similarity in JoStARS
#' @param method how should the similarity between the prec matrices be? Symmetric by default, meaning all are equally different. If not, one will stand out as completely unrelated.
#' @param include.jostars should we perform JoStARS?
#' @param include.JGL should the joint graphical lasso be included?
#' @param stars.thresh the variability threshold to use in the graphical lasso tuned by StARS
#' @param var.thesh.stars if include.jostars, the variability threshold to use. Default \eqn{0.1}
#' @param penalize.diagonal should the diagonal be penalized in the graphical lasso-based methods?
#' @param verbose logical indicator of printing information at each iteration
#' @param scale should the data be scaled?
#' @param save.res.jointGHS should Lambda, Nu and Theta be saved for all N simulations? Requires a lot of memory. Default FALSE.
#' @return simulation results, including sparsity, precision, recall and specificity
perform_jointGHS_simulation = function(K, n.vals, p, N=100, seeds=sample(1:1000,N), nCores = 3, frac.disagreement = 0, tau_sq = 10, tau_sq_ghs = 1, ebic.gamma = 0.2, method='symmetric', 
                                       include.jostars=TRUE, include.JGL = FALSE, stars.thresh=0.05, var.thresh.jostars = 0.05,penalize.diagonal=FALSE ,verbose=TRUE, scale=TRUE, 
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
  
  # JGL results (tuned by AIC)
  res$opt.sparsities.jgl = matrix(0,N,K)
  res$precisions.jgl =  matrix(0,N,K)
  res$specificities.jgl =  matrix(0,N,K)
  res$recalls.jgl =  matrix(0,N,K)
  res$matrix.distances.jgl =  matrix(0,N,K)
  
  # JoStARS results
  res$opt.sparsities.jostars = matrix(0,N,K)
  res$precisions.jostars =  matrix(0,N,K)
  res$specificities.jostars =  matrix(0,N,K)
  res$recalls.jostars =  matrix(0,N,K)
  res$matrix.distances.jostars =  matrix(0,N,K)
  
  # GHS results
  res$opt.sparsities.ghs = matrix(0,N,K)
  res$precisions.ghs =  matrix(0,N,K)
  res$specificities.ghs =  matrix(0,N,K)
  res$recalls.ghs =  matrix(0,N,K)
  res$matrix.distances.ghs =  matrix(0,N,K)
  
  # Glasso results (tuned by AIC)
  res$opt.sparsities.glasso = matrix(0,N,K)
  res$precisions.glasso =  matrix(0,N,K)
  res$specificities.glasso =  matrix(0,N,K)
  res$recalls.glasso =  matrix(0,N,K)
  res$matrix.distances.glasso =  matrix(0,N,K)
  
  # Start by generating the precision matrices 
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
  if(scale) cov.matrices[[1]] = cov2cor(cov.matrices[[1]])
  if(scale) prec.matrices[[1]] = cov2cor(prec.matrices[[1]])
  
  # Avoid rounding errors leading to matrices not being symmetric
  if(!matrixcalc::is.symmetric.matrix(cov.matrices[[1]])){
    cov.matrices[[1]] = round(cov.matrices[[1]],8)
  }
  if(method=='symmetric'){
    for(k in 2:K){
      huge.tmp = mutate.graph(huge.init,frac.disagreement,scale)
      cov.matrices[[k]] = huge.tmp$cov.mat
      prec.matrices[[k]] = huge.tmp$prec.mat
      # Avoid rounding errors leading to matrices not being symmetric
      if(!matrixcalc::is.symmetric.matrix(cov.matrices[[k]])){
        cov.matrices[[k]] = round(cov.matrices[[k]],8)
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
                                     ebic.gamma=ebic.gamma,tau_sq = tau_sq, tau_sq_ghs=tau_sq_ghs,include.jostars=include.jostars, include.JGL=include.JGL, 
                                     stars.thresh=stars.thresh, var.thresh.jostars=var.thresh.jostars, penalize.diagonal=penalize.diagonal,seed=seeds[i], 
                                     save.res.jointGHS = save.res.jointGHS);
  }
  registerDoSEQ()
  
  # Save results from each replicate
  for(i in 1:N){
    est.tmp = res.list[[i]]
    # Results from jointGHS
    res$opt.sparsities[i,] = est.tmp$opt.sparsities 
    res$matrix.distances[i,] = est.tmp$matrix.distances
    res$precisions[i,] = est.tmp$precisions
    res$recalls[i,] = est.tmp$recalls
    res$specificities[i,] =  est.tmp$specificities
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
    # Results from jostars
    if(include.jostars){
      res$opt.sparsities.jostars[i,] = est.tmp$opt.sparsities.jostars
      res$matrix.distances.jostars[i,] = est.tmp$matrix.distances.jostars
      res$precisions.jostars[i,] = est.tmp$precisions.jostars
      res$recalls.jostars[i,] = est.tmp$recalls.jostars
      res$specificities.jostars[i,] =  est.tmp$specificities.jostars
    }
    # Results from ghs
    res$opt.sparsities.ghs[i,] = est.tmp$opt.sparsities.ghs
    res$matrix.distances.ghs[i,] = est.tmp$matrix.distances.ghs
    res$precisions.ghs[i,] = est.tmp$precisions.ghs
    res$recalls.ghs[i,] = est.tmp$recalls.ghs
    res$specificities.ghs[i,] =  est.tmp$specificities.ghs
    # Results from glasso
    res$matrix.distances.glasso[i,] = est.tmp$matrix.distances.glasso
    res$precisions.glasso[i,] = est.tmp$precisions.glasso
    res$recalls.glasso[i,] = est.tmp$recalls.glasso
    res$specificities.glasso[i,] =  est.tmp$specificities.glasso
    res$opt.sparsities.glasso[i,] = est.tmp$opt.sparsities.glasso
  }
  # Mean results from jointGHS
  res$mean.opt.sparsities= colMeans(res$opt.sparsities)
  res$mean.precisions= colMeans(res$precisions)
  res$mean.recalls = colMeans(res$recalls)
  res$mean.specificities =  colMeans(res$specificities)
  res$mean.matrix.distances = colMeans(res$matrix.distances)
  
  # Mean results from JGL
  if(include.JGL){
    res$mean.opt.sparsities.jgl = colMeans(res$opt.sparsities.jgl)
    res$mean.precisions.jgl = colMeans(res$precisions.jgl)
    res$mean.recalls.jgl = colMeans(res$recalls.jgl)
    res$mean.specificities.jgl =  colMeans(res$specificities.jgl)
    res$mean.matrix.distances.jgl = colMeans(res$matrix.distances.jgl)
  }
  
  # Mean results from JoStARS
  if(include.jostars){
    res$mean.opt.sparsities.jostars = colMeans(res$opt.sparsities.jostars)
    res$mean.precisions.jostars = colMeans(res$precisions.jostars)
    res$mean.recalls.jostars = colMeans(res$recalls.jostars)
    res$mean.specificities.jostars =  colMeans(res$specificities.jostars)
    res$mean.matrix.distances.jostars = colMeans(res$matrix.distances.jostars)
  }
  
  # Mean results from GHS
  res$mean.opt.sparsities.ghs = colMeans(res$opt.sparsities.ghs)
  res$mean.precisions.ghs = colMeans(res$precisions.ghs)
  res$mean.recalls.ghs = colMeans(res$recalls.ghs)
  res$mean.specificities.ghs =  colMeans(res$specificities.ghs)
  res$mean.matrix.distances.ghs = colMeans(res$matrix.distances.ghs)
  
  # Mean results from glasso
  res$mean.opt.sparsities.glasso= colMeans(res$opt.sparsities.glasso)
  res$mean.precisions.glasso= colMeans(res$precisions.glasso)
  res$mean.recalls.glasso = colMeans(res$recalls.glasso)
  res$mean.specificities.glasso =  colMeans(res$specificities.glasso)
  res$mean.matrix.distances.glasso = colMeans(res$matrix.distances.glasso)
  
  res$true.prec.matrices = prec.matrices
  res$true.sparsity = spars.init
  return(res)
  
}

# Function for performing one iteration -----------------------------------------

jointGHS_simulation_one_iteration = function(n.vals,cov.matrices,prec.matrices,scale,ebic.gamma,tau_sq, tau_sq_ghs,include.jostars,include.JGL, 
                                             stars.thresh, var.thresh.jostars, penalize.diagonal,seed, save.res.jointGHS) {
  y = list()
  K=length(n.vals)
  p=ncol(prec.matrices[[1]])
  glasso.res=list()
  ghs.res = list()
  set.seed(seed)
  # Generate data. 
  for(k in 1:K){
    y[[k]] = mvtnorm::rmvnorm(n.vals[k], mean=rep(0,p), cov.matrices[[k]])
    if (scale) y[[k]] = scale(y[[k]])
    # Use single-network methods
    glasso.tmp = huge(y[[k]],method='glasso',verbose = F)
    glasso.res[[k]] = huge.select(glasso.tmp,criterion='stars', stars.thresh = stars.thresh, verbose = F)
    ghs.res[[k]] = fastGHS(X=y[[k]], tau_sq = tau_sq_ghs, verbose=FALSE, fix_tau = TRUE)$theta
    ghs.res[[k]][which(abs(ghs.res[[k]])<1e-5, arr.ind=T)] = 0
  }
  # Perform joint methods
  res.full.tmp = jointGHS(X=y, tau_sq = tau_sq, verbose = F, epsilon=1e-3,fix_tau = T)
  res.tmp = res.full.tmp$theta
  
  if(include.JGL){
    jgl.tmp = JGL_select_AIC(Y=y,penalty='fused',nlambda1=10,lambda1.min=0.01,lambda1.max=1,nlambda2=10,lambda2.min=0,lambda2.max=0.1,lambda2.init=0.01,
                             penalize.diagonal=penalize.diagonal)
  }
  if(include.jostars){
    jostars.tmp = JoStARS::JoStARS(Y=y,var.thresh = var.thresh.jostars,scale=T,nlambda1=20,
                                   lambda1.min=0.01,lambda1.max=1, nlambda2=20,lambda2.min=0,lambda2.max = 0.1, lambda2.init = 0.01,
                                   ebic.gamma=ebic.gamma,verbose=F,parallelize = FALSE, penalize.diagonal=penalize.diagonal)
  }
  est=list()
  
  # If save.res.jointGHS==T, save Lambda, expectation of 1/Nu and Theta
  if(save.res.jointGHS){
    est$E_NuInv = res.full.tmp$E_NuInv # A K by K matrix
    est$theta = res.full.tmp$theta # A list of length K
    est$Lambda_sq = res.full.tmp$Lambda_sq # A list of length K
  }
  
  # Results from jointGHS
  res.tmp = lapply(res.tmp, cov2cor)
  est$opt.sparsities  = unlist(lapply(res.tmp, FUN= function(m) tailoredGlasso::sparsity(abs(res.tmp[[k]])>1e-5)))
  est$matrix.distances = sapply(1:K,FUN=function(k) matrix.distance(prec.matrices[[k]], res.tmp[[k]]))
  est$precisions =  sapply(1:K,FUN=function(k) precision(prec.matrices[[k]]!=0, abs(res.tmp[[k]])>1e-5 ))
  est$recalls =  sapply(1:K,FUN=function(k) recall(prec.matrices[[k]]!=0, abs(res.tmp[[k]])>1e-5 ))
  est$specificities =  sapply(1:K,FUN=function(k) specificity(prec.matrices[[k]]!=0, abs(res.tmp[[k]])>1e-5 ))
  
  # Results from JGL
  if(include.JGL){
    est$opt.sparsities.jgl = jgl.tmp$opt.sparsities
    est$matrix.distances.jgl = sapply(1:K,FUN=function(k) matrix.distance(prec.matrices[[k]], jgl.tmp$opt.fit[[k]]))
    est$precisions.jgl = sapply(1:K,FUN=function(k) precision(prec.matrices[[k]]!=0, jgl.tmp$opt.fit[[k]]!=0))
    est$recalls.jgl = sapply(1:K,FUN=function(k) recall(prec.matrices[[k]]!=0, jgl.tmp$opt.fit[[k]]!=0))
    est$specificities.jgl = sapply(1:K,FUN=function(k) specificity(prec.matrices[[k]]!=0, jgl.tmp$opt.fit[[k]]!=0))
  }
  
  # Results from JoStARS
  if(include.jostars){
    est$opt.sparsities.jostars = jostars.tmp$opt.sparsities
    est$matrix.distances.jostars = sapply(1:K,FUN=function(k) matrix.distance(prec.matrices[[k]], jostars.tmp$opt.fit[[k]]))
    est$precisions.jostars = sapply(1:K,FUN=function(k) precision(prec.matrices[[k]]!=0, jostars.tmp$opt.fit[[k]]!=0))
    est$recalls.jostars = sapply(1:K,FUN=function(k) recall(prec.matrices[[k]]!=0, jostars.tmp$opt.fit[[k]]!=0))
    est$specificities.jostars = sapply(1:K,FUN=function(k) specificity(prec.matrices[[k]]!=0, jostars.tmp$opt.fit[[k]]!=0))
  }
  
  # Results from GHS
  est$opt.sparsities.ghs = unlist(lapply(ghs.res, FUN= function(m) tailoredGlasso::sparsity(m!=0)))
  est$matrix.distances.ghs = sapply(1:K,FUN=function(k) matrix.distance(prec.matrices[[k]], ghs.res[[k]]))
  est$precisions.ghs = sapply(1:K,FUN=function(k) precision(prec.matrices[[k]]!=0, ghs.res[[k]]!=0))
  est$recalls.ghs = sapply(1:K,FUN=function(k) recall(prec.matrices[[k]]!=0, ghs.res[[k]]!=0))
  est$specificities.ghs =  sapply(1:K,FUN=function(k) specificity(prec.matrices[[k]]!=0, ghs.res[[k]]!=0))
  
  # Results from glasso
  est$matrix.distances.glasso = sapply(1:K,FUN=function(k) matrix.distance(prec.matrices[[k]], glasso.res[[k]]$opt.icov))
  est$precisions.glasso = sapply(1:K,FUN=function(k) precision(prec.matrices[[k]]!=0, glasso.res[[k]]$opt.icov!=0))
  est$recalls.glasso = sapply(1:K,FUN=function(k) recall(prec.matrices[[k]]!=0, glasso.res[[k]]$opt.icov!=0))
  est$specificities.glasso =  sapply(1:K,FUN=function(k) specificity(prec.matrices[[k]]!=0, glasso.res[[k]]$opt.icov!=0))
  est$opt.sparsities.glasso = unlist(lapply(glasso.res,FUN=function(l) l$opt.sparsity))
  return(est)
}




