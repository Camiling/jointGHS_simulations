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
source('GHS/GHS.R')

#' Perform jointGHS simulations
#' 
#' This function performs simulations for jointGHS, averaging over the results. 
#' 
#' @param n the number of observations in the data set
#' @param p the number of nodes
#' @param N the number of simulations to perform
#' @param seeds seeds to use in each of the \eqn{N} simulations. A vector of length \eqn{N}. 
#' @param nCores how many cores should be used
#' @param ebic.gamma the additional penalty term in the extended BIC criterion for selecting similarity in stabJGL
#' @param include.glasso should we perform  the graphical lasso?
#' @param include.GHS should the Gibbs GHS be performed?
#' @param include.fastGHS should the ECMGHS be performed?
#' @param stars.thresh the variability threshold to use in the graphical lasso tuned by StARS
#' @param penalize.diagonal should the diagonal be penalized in the graphical lasso-based methods?
#' @param verbose logical indicator of printing information at each iteration
#' @param scale should the data be scaled?
#' @return simulation results, including sparsity, precision, recall and specificity
perform_fastGHS_simulation = function(n, p, N=100, seeds=sample(1:1000,N), nCores = 3, ebic.gamma = 0.2,
                                       include.glasso = TRUE, include.GHS=TRUE, include.fastGHS=TRUE,
                                       stars.thresh=0.03,penalize.diagonal=FALSE ,verbose=TRUE, scale=TRUE, theta.init=T,save_Q=F){
  
  res=list()
  # Glasso results
  if(include.glasso){
  res$opt.sparsities.glasso = rep(0,N)
  res$precisions.glasso =  rep(0,N)
  res$specificities.glasso =  rep(0,N)
  res$recalls.glasso =  rep(0,N)
  res$matrix.distances.glasso = rep(0,N)
  }
  # GHS results
  if(include.GHS){
    res$opt.sparsities.ghs = rep(0,N)
    res$precisions.ghs =  rep(0,N)
    res$specificities.ghs =  rep(0,N)
    res$recalls.ghs =  rep(0,N)
    res$matrix.distances.ghs = rep(0,N)
  }
  # fastGHS results 
  if(include.fastGHS){
    res$opt.sparsities.fastghs = rep(0,N)
    res$precisions.fastghs =  rep(0,N)
    res$specificities.fastghs =  rep(0,N)
    res$recalls.fastghs =  rep(0,N)
    res$matrix.distances.fastghs = rep(0,N)
    if(save_Q){
      res$Q_vals = list()
    }
  }
  
  # Generate the precision matrix
  huge.init = huge.generator(n,p,graph='scale-free',verbose = F,v=0.5,u=0.05)
  theta.init = huge.init$omega
  theta.init[which(abs(theta.init)<1e-5,arr.ind=T)] = 0
  spars.init = huge.init$sparsity
  cov.mat = huge.init$sigma
  prec.mat = theta.init
  
  # Scale?
  if(scale) cov.mat = cov2cor(cov.mat)
  if(scale) prec.mat = cov2cor(prec.mat)
  
  # Avoid rounding errors leading to matrices not being symmetric
  if(!matrixcalc::is.symmetric.matrix(cov.mat)){
    cov.mat = round(cov.mat,8)
  }

  registerDoParallel(nCores)
  res.list = foreach (i=1:N) %dopar% {
    jointGHS_simulation_one_iteration(n=n,cov.mat=cov.mat,prec.mat=prec.mat,scale=scale,
                                      ebic.gamma=ebic.gamma,include.glasso=include.glasso, include.GHS=include.GHS, include.fastGHS=include.fastGHS, 
                                      stars.thresh=stars.thresh,penalize.diagonal=penalize.diagonal,seed=seeds[i],theta.init=theta.init,save_Q=save_Q);
  }
  registerDoSEQ()
  # Save results from each replicate
  for(i in 1:N){
    est.tmp = res.list[[i]]
    # Results from GHS
    if(include.GHS){
      res$opt.sparsities.ghs[i] = est.tmp$opt.sparsities.ghs
      res$matrix.distances.ghs[i] = est.tmp$matrix.distances.ghs
      res$precisions.ghs[i] = est.tmp$precisions.ghs
      res$recalls.ghs[i] = est.tmp$recalls.ghs
      res$specificities.ghs[i] =  est.tmp$specificities.ghs
    }
    # Results from fastGHS
    if(include.fastGHS){
      res$opt.sparsities.fastghs[i] = est.tmp$opt.sparsities.fastghs
      res$matrix.distances.fastghs[i] = est.tmp$matrix.distances.fastghs
      res$precisions.fastghs[i] = est.tmp$precisions.fastghs
      res$recalls.fastghs[i] = est.tmp$recalls.fastghs
      res$specificities.fastghs[i] =  est.tmp$specificities.fastghs
      if(save_Q){
        res$Q_vals[[i]] = est.tmp$Q_vals 
      }
    }
    # Results from glasso
    if(include.glasso){
      res$matrix.distances.glasso[i] = est.tmp$matrix.distances.glasso
      res$precisions.glasso[i] = est.tmp$precisions.glasso
      res$recalls.glasso[i] = est.tmp$recalls.glasso
      res$specificities.glasso[i] =  est.tmp$specificities.glasso
      res$opt.sparsities.glasso[i] = est.tmp$opt.sparsities.glasso      
    }
  }
  # Mean results from GHS
  if(include.GHS){
    res$mean.opt.sparsities.ghs = mean(res$opt.sparsities.ghs)
    res$mean.precisions.ghs = mean(res$precisions.ghs)
    res$mean.recalls.ghs = mean(res$recalls.ghs)
    res$mean.specificities.ghs =  mean(res$specificities.ghs)
    res$mean.matrix.distances.ghs = mean(res$matrix.distances.ghs)    
  }
  # Mean results from fastGHS
  if(include.fastGHS){
    res$mean.opt.sparsities.fastghs = mean(res$opt.sparsities.fastghs)
    res$mean.precisions.fastghs = mean(res$precisions.fastghs)
    res$mean.recalls.fastghs = mean(res$recalls.fastghs)
    res$mean.specificities.fastghs =  mean(res$specificities.fastghs)
    res$mean.matrix.distances.fastghs = mean(res$matrix.distances.fastghs)    
  }
  # Mean results from glasso
  if(include.glasso){
    res$mean.opt.sparsities.glasso= mean(res$opt.sparsities.glasso)
    res$mean.precisions.glasso= mean(res$precisions.glasso)
    res$mean.recalls.glasso = mean(res$recalls.glasso)
    res$mean.specificities.glasso =  mean(res$specificities.glasso)
    res$mean.matrix.distances.glasso = mean(res$matrix.distances.glasso)    
  }
  res$true.prec.mat = prec.mat
  res$true.sparsity = spars.init
  return(res)
  
}

# Function for performing one iteration -----------------------------------------

jointGHS_simulation_one_iteration = function(n,cov.mat,prec.mat,scale,ebic.gamma,include.glasso,include.GHS, include.fastGHS, 
                                             stars.thresh, penalize.diagonal,seed,theta.init=F,save_Q=F) {
  p=ncol(prec.mat)
  set.seed(seed)
  # Generate data. 
  y = mvtnorm::rmvnorm(n, mean=rep(0,p), cov.mat)
  if (scale) y= scale(y)
  # Infer networks
  if(include.glasso){
  glasso.tmp = huge(y,method='glasso',verbose = F)
  glasso.res = huge.select(glasso.tmp,criterion='stars', stars.thresh = stars.thresh, verbose = F)    
  }
  if(include.fastGHS){
    if(theta.init){
      huge.init.random = huge.generator(n,p,graph='scale-free',verbose = F,v=0.5,u=0.05)
      theta.init.random = cov2cor(huge.init.random$omega)
      if(!matrixcalc::is.symmetric.matrix(theta.init.random)){
        theta.init.random = round(theta.init.random,8)
      }
      fastghs.res = fastGHS::fastGHS(X=y, AIC_selection = T, AIC_eps = 0.1, epsilon = 1e-3, verbose=FALSE, theta=theta.init.random,save_Q=save_Q)
      Q_vals = c(fastghs.res$Q_vals)
      fastghs.res = fastghs.res$theta
      fastghs.res[which(abs(fastghs.res)<1e-5, arr.ind=T)] = 0
    }
    else{
      fastghs.res = fastGHS::fastGHS(X=y, AIC_selection = T, AIC_eps = 0.1, epsilon = 1e-3, verbose=FALSE)$theta
      fastghs.res[which(abs(fastghs.res)<1e-5, arr.ind=T)] = 0     
    }

  }
  if(include.GHS){
    invalid=T
    while(invalid){
      ghs.res = tryCatch(GHS(t(y)%*%y,n,burnin=100,nmc=1000), error=function(s) NULL)
      if(is.null(ghs.res)) {
        invalid=T
      }
      else{
        theta.est.ghs = cov2cor(apply(ghs.res$thetas.sampled, c(1,2), mean))
        # Threshold to get comparable sparsity
        theta.est.off.diag.ghs = theta.est.ghs
        diag(theta.est.off.diag.ghs) = NA
        theta.est.ghs[which(abs(theta.est.ghs) < quantile(abs(theta.est.off.diag.ghs),1-tailoredGlasso::sparsity(fastghs.res!=0),na.rm = T), arr.ind = T)] = 0
        ghs.res = theta.est.ghs
        invalid=F
      }
    }
  }
  est=list()
  # Results from GHS
  if(include.GHS){
    est$opt.sparsities.ghs = tailoredGlasso::sparsity(ghs.res!=0)
    est$matrix.distances.ghs = matrix.distance(prec.mat, ghs.res)
    est$precisions.ghs = precision(prec.mat!=0, ghs.res!=0)
    est$recalls.ghs = recall(prec.mat!=0, ghs.res!=0)
    est$specificities.ghs = specificity(prec.mat!=0, ghs.res!=0)
  }
  # Results from fastGHS
  if(include.fastGHS){
    est$opt.sparsities.fastghs = tailoredGlasso::sparsity(fastghs.res!=0)
    est$matrix.distances.fastghs = matrix.distance(prec.mat, fastghs.res)
    est$precisions.fastghs = precision(prec.mat!=0, fastghs.res!=0)
    est$recalls.fastghs = recall(prec.mat!=0, fastghs.res!=0)
    est$specificities.fastghs =  specificity(prec.mat!=0, fastghs.res!=0)
    if(save_Q){
      est$Q_vals = Q_vals
    }
  }
  # Results from glasso
  if(include.glasso){
    est$matrix.distances.glasso = matrix.distance(prec.mat, glasso.res$opt.icov)
    est$precisions.glasso = precision(prec.mat!=0, glasso.res$opt.icov!=0)
    est$recalls.glasso =  recall(prec.mat!=0, glasso.res$opt.icov!=0)
    est$specificities.glasso =  specificity(prec.mat!=0, glasso.res$opt.icov!=0)
    est$opt.sparsities.glasso = glasso.res$opt.sparsity
  }
  return(est)
}




