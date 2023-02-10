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
source('GemBag/GemBag_implementation/BIC_GemBag.R')
Rcpp::sourceCpp('GemBag/GemBag_implementation/GemBag-algo.cpp')

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
#' @param verbose logical indicator of printing information at each iteration
#' @param scale should the data be scaled?
#' @param save.edgeagreement should the edge agreement of the estimates be saved?
#' @return simulation results, including sparsity, precision, recall and specificity
perform_jointGHS_simulation_gembag = function(K, n.vals, p, N=100, seeds=sample(1:1000,N), nCores = 3, frac.disagreement = 0, method='symmetric', verbose=TRUE, scale=TRUE, 
                                              save.edgeagreement=FALSE){
  
  res=list()
  # GemBag results 
  res$opt.sparsities.gembag = matrix(0,N,K)
  res$precisions.gembag =  matrix(0,N,K)
  res$specificities.gembag =  matrix(0,N,K)
  res$recalls.gembag =  matrix(0,N,K)
  res$matrix.distances.gembag =  matrix(0,N,K)
  res$edge.agreement.gembag = rep(0,K)
  res$n.edges.est.gembag =  matrix(0,N,K)
  
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
        huge.tmp = mutate.graph(huge.init,frac.disagreement,scale, generate.data = F, larger.partialcor = T)
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
    jointGHS_simulation_one_iteration_gembag(n.vals=n.vals,cov.matrices=cov.matrices,prec.matrices=prec.matrices,scale=scale,seed=seeds[i], save.edgeagreement=save.edgeagreement);
  }
  registerDoSEQ()
  
  # Save results from each replicate
  for(i in 1:N){
    est.tmp = res.list[[i]]
    # Results from gembag
    res$opt.sparsities.gembag[i,] = est.tmp$opt.sparsities.gembag
    res$matrix.distances.gembag[i,] = est.tmp$matrix.distances.gembag
    res$precisions.gembag[i,] = est.tmp$precisions.gembag
    res$recalls.gembag[i,] = est.tmp$recalls.gembag
    res$specificities.gembag[i,] =  est.tmp$specificities.gembag 
    if(save.edgeagreement){
      res$edge.agreement.gembag[i] = est.tmp$edge.agreement.gembag
      res$n.edges.est.gembag[i,] =  est.tmp$n.edges.est.gembag
    }
  }

  # Mean results from gembag
  res$mean.opt.sparsities.gembag = colMeans(res$opt.sparsities.gembag)
  res$mean.precisions.gembag = colMeans(res$precisions.gembag)
  res$mean.recalls.gembag = colMeans(res$recalls.gembag)
  res$mean.specificities.gembag =  colMeans(res$specificities.gembag)
  res$mean.matrix.distances.gembag = colMeans(res$matrix.distances.gembag)
  res$true.prec.matrices = prec.matrices
  res$true.sparsity = spars.init
  if(save.edgeagreement){
    res$mean.edge.agreement.gembag = mean(res$edge.agreement.gembag)
    res$mean.n.edges.est.gembag =  colMeans(res$n.edges.est.gembag)
  }
  return(res)
  
}

# Function for performing one iteration -----------------------------------------

jointGHS_simulation_one_iteration_gembag = function(n.vals,cov.matrices,prec.matrices,scale,seed, save.edgeagreement) {
  y = list()
  K=length(n.vals)
  p=ncol(prec.matrices[[1]])
  set.seed(seed)
  S_l = list()
  # Generate data. 
  for(k in 1:K){
    y[[k]] = mvtnorm::rmvnorm(n.vals[k], mean=rep(0,p), cov.matrices[[k]])
    if (scale) y[[k]] = scale(y[[k]])
    S_l[[k]] = cov(y[[k]])
  }
  # Perform joint method
  # Set of hyperparameters
  v0_l <- c(0.25, 0.5, 0.75, 1) * sqrt(1/n.vals[1]/log(p))
  v1_l <- c(2.5, 5, 7.5, 10) * sqrt(1/n.vals[1]/log(p))
  p1 <- 0.5
  # Tuning by BIC
  hyper <- Tune_GemBag(v0_l, v1_l,S_l, n.vals, maxiter=20, p1, p_2=0.5)
  v0 <- hyper$v0
  v1 <- hyper$v1
  # Final estimate
  gembag.tmp <- GemBag(S_l=S_l, n= n.vals, 
                  v_0=v0, v_1=v1, tau=v0, 
                  p_1=p1, p_2=1,
                  maxiter=20)
  names(gembag.tmp) <- c('Theta', 'P', 'W')
  gembag.tmp = gembag.tmp$Theta
  est=list()
  
  # Results from gembag
  gembag.tmp = lapply(gembag.tmp, cov2cor)
  est$opt.sparsities.gembag  = unlist(lapply(gembag.tmp, FUN= function(m) tailoredGlasso::sparsity(m!=0)))
  est$matrix.distances.gembag = sapply(1:K,FUN=function(k) matrix.distance(prec.matrices[[k]], gembag.tmp[[k]]))
  est$precisions.gembag = sapply(1:K,FUN=function(k) precision(prec.matrices[[k]]!=0, gembag.tmp[[k]]!=0))
  est$recalls.gembag = sapply(1:K,FUN=function(k) recall(prec.matrices[[k]]!=0, gembag.tmp[[k]]!=0))
  est$specificities.gembag = sapply(1:K,FUN=function(k) specificity(prec.matrices[[k]]!=0, gembag.tmp[[k]]!=0))
  if(save.edgeagreement){
    est$edge.agreement.gembag = (sum(apply(simplify2array(gembag.tmp), c(1,2), FUN = function(m) sum(m!=0)==K))-p)/2 # How many edges common to all K nets?
    est$n.edges.est.gembag =  sapply(1:K,FUN=function(k) (sum(gembag.tmp[[k]]!=0)-p)/2)
  }
  return(est)
}

print_results_jointGHS_show_SD_gembag = function(obj.list,fracs.mutated, show.distance=F,show.specificity=F, show.edgedisagreement=F){
  # obj is a list of objects returned by perform_jointGHS_simulation.
  # fracs.mutated is a vector of the mutated fraction in each simulation object
  # show.distance: should the matrix distance be printed?
  # show.specificity: should the matrix distance be printed?
  # Function for printing mean sparsities, precisions, recalls and matrix distances when several data sets were generated.
  # Note that we print the results for the different graphs on the same lines. 
  K = length(obj.list[[1]]$mean.opt.sparsities)
  # Loop over each scenario
  for (i in 1:length(obj.list)){
    obj=obj.list[[i]]
    cat(fracs.mutated[i])
    cat('  & GemBag ')  
    if(show.edgedisagreement){
      cat('&', round((1-2*obj$mean.edge.agreement.gembag/ (sum(obj$mean.n.edges.est.gembag)))*100))
    }
      for(k in 1:K){
        cat(' && ',round(obj$mean.opt.sparsities.gembag[k],3), '(',round(sd(obj$opt.sparsities.gembag[,k]),3),')',' & ',
        round(obj$mean.precisions.gembag[k],2),'(',round(sd(obj$precisions.gembag[,k]),2),')',' & ',
        round(obj$mean.recalls.gembag[k],2),'(',round(sd(obj$recalls.gembag[,k]),2),')')
      if(show.specificity)cat('&',round(obj$mean.specificities.gembag[k],2), '(',round(sd(obj$specificities.gembag[,k]),2),')') 
      if(show.distance) cat(' & ',round(obj$mean.matrix.distances.gembag[k],3))
      }
    cat(' \\\\ \n \\hline \n')
  }
}


