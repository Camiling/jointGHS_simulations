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
source('GemBag/Gembag_implementation/BIC_GemBag.R')
Rcpp::sourceCpp('GemBag/Gembag_implementation/GemBag-algo.cpp')

#' Perform jointGHS simulations, for creating ROC curves
#' 
#' This function performs simulations for joint methods, averaging over the results. 
#' 
#' @param K number of data sets/networks
#' @param n.vals the number of observations in each data set. A vector
#' @param p the number of nodes
#' @param N the number of simulations to perform
#' @param seeds seeds to use in each of the \eqn{N} simulations. A vector of length \eqn{N}. 
#' @param nCores how many cores should be used
#' @param frac.disagreement the fraction of edges that the networks should disagree on
#' @param method Which joint method to use?
#' @param penalize.diagonal should the diagonal be penalized in the graphical lasso-based methods?
#' @param verbose logical indicator of printing information at each iteration
#' @param scale should the data be scaled?
#' @param save.edgeagreement should the number of edges the estimates agree on be saved?
#' @return simulation results, including sparsity, precision, recall and specificity
perform_jointGHS_simulation_ROC = function(K, n.vals, p, N=100, method = 'jointGHS', seeds=sample(1:1000,N), nCores = 3, frac.disagreement = 0, 
                                           penalize.diagonal=FALSE ,verbose=TRUE, scale=TRUE, save.edgeagreement=F){
  
  res=list()
  include.jointGHS=FALSE
  include.JGL= FALSE
  include.SSJGL=FALSE
  include.gembag=FALSE
  if(method=='jointGHS'){
    include.jointGHS=TRUE
  }
  else if(method=='JGL'){
    include.JGL=TRUE
  }
  else if(method=='SSJGL'){
    include.SSJGL=TRUE
  }
  else if(method=='GemBag'){
    include.gembag=TRUE
  }
  # JointGHS results
  n.tau.jointGHS = 500
  res$precisions =  array(0,c(N,K,n.tau.jointGHS))
  res$recalls =  array(0,c(N,K,n.tau.jointGHS))
  res$TPR =  array(0,c(N,K,n.tau.jointGHS))
  res$FPR=  array(0,c(N,K,n.tau.jointGHS))
  
  # SSJGL results 
  n.lam.eff.ssjgl = length(seq(2,6,by=0.1))
  res$precisions.ssjgl =  array(0,c(N,K,n.lam.eff.ssjgl))
  res$recalls.ssjgl =  array(0,c(N,K,n.lam.eff.ssjgl))
  res$TPR.ssjgl =  array(0,c(N,K,n.lam.eff.ssjgl))
  res$FPR.ssjgl =  array(0,c(N,K,n.lam.eff.ssjgl))
  
  # JGL results 
  n.lam1.jgl = 50
  res$precisions.jgl =  array(0,c(N,K,n.lam1.jgl))
  res$recalls.jgl =  array(0,c(N,K,n.lam1.jgl))
  res$TPR.jgl =  array(0,c(N,K,n.lam1.jgl))
  res$FPR.jgl =  array(0,c(N,K,n.lam1.jgl))
  
  # GemBag results 
  n.thresh.gembag = 1000
  res$precisions.gembag =  array(0,c(N,K,n.thresh.gembag))
  res$recalls.gembag =  array(0,c(N,K,n.thresh.gembag))
  res$TPR.gembag =  array(0,c(N,K,n.thresh.gembag))
  res$FPR.gembag =  array(0,c(N,K,n.thresh.gembag))
  
  
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
  registerDoParallel(nCores)
  res.list = foreach (i=1:N) %dopar% {
    jointGHS_simulation_one_iteration_ROC(n.vals=n.vals,cov.matrices=cov.matrices,prec.matrices=prec.matrices,scale=scale,
                                          include.jointGHS=include.jointGHS,include.SSJGL=include.SSJGL,include.JGL=include.JGL,
                                          include.gembag = include.gembag, penalize.diagonal=penalize.diagonal,seed=seeds[i], n.tau.jointGHS=n.tau.jointGHS, 
                                          n.lam.eff.ssjgl=n.lam.eff.ssjgl,n.lam1.jgl=n.lam1.jgl,n.thresh.gembag=n.thresh.gembag);
  }
  registerDoSEQ()
  
  # Save results from each replicate
  for(i in 1:N){
    est.tmp = res.list[[i]]
    
    # Results from jointGHS
    if(include.jointGHS){
      res$precisions[i,,] = est.tmp$precisions
      res$recalls[i,,] = est.tmp$recalls
      res$TPR[i,,] = est.tmp$TPR
      res$FPR[i,,] = est.tmp$FPR
    }
    # Results from jgl
    else if(include.JGL){
      res$precisions.jgl[i,,] = est.tmp$precisions.jgl
      res$recalls.jgl[i,,] = est.tmp$recalls.jgl
      res$TPR.jgl[i,,] = est.tmp$TPR.jgl
      res$FPR.jgl[i,,] = est.tmp$FPR.jgl
    }
    # Results from ssjgl
    else if(include.SSJGL){
      res$precisions.ssjgl[i,,] = est.tmp$precisions.ssjgl
      res$recalls.ssjgl[i,,] = est.tmp$recalls.ssjgl
      res$TPR.ssjgl[i,,] = est.tmp$TPR.ssjgl
      res$FPR.ssjgl[i,,] = est.tmp$FPR.ssjgl
    }
    # Results from GemBag
    else if(include.gembag){
      res$precisions.gembag[i,,] = est.tmp$precisions.gembag
      res$recalls.gembag[i,,] = est.tmp$recalls.gembag
      res$TPR.gembag[i,,] = est.tmp$TPR.gembag
      res$FPR.gembag[i,,] = est.tmp$FPR.gembag
    }
  }
  # Mean results from jointGHS
  if(include.jointGHS){
    res$mean.precisions=apply(res$precisions, 3, colMeans) # Becomes K x n.tau.jointGHS
    res$mean.recalls = apply(res$recalls, 3, colMeans) 
    res$mean.TPR =apply(res$TPR, 3, colMeans) 
    res$mean.FPR = apply(res$FPR, 3, colMeans) 
  }
  # Mean results from JGL
  if(include.JGL){
    res$mean.precisions.jgl = apply(res$precisions.jgl, 3, colMeans)
    res$mean.recalls.jgl = apply(res$recalls.jgl, 3, colMeans)
    res$mean.TPR.jgl =apply(res$TPR.jgl, 3, colMeans) 
    res$mean.FPR.jgl = apply(res$FPR.jgl, 3, colMeans) 
  }
  # Mean results from SSJGL
  if(include.SSJGL){
    res$mean.precisions.ssjgl = apply(res$precisions.ssjgl, 3, colMeans) 
    res$mean.recalls.ssjgl = apply(res$recalls.ssjgl, 3, colMeans) 
    res$mean.TPR.ssjgl =apply(res$TPR.ssjgl, 3, colMeans) 
    res$mean.FPR.ssjgl = apply(res$FPR.ssjgl, 3, colMeans) 
  }
  # Mean results from GemBag
  if(include.gembag){
    res$mean.precisions.gembag = apply(res$precisions.gembag, 3, colMeans)
    res$mean.precisions.gembag = cbind(rep(1,K),res$mean.precisions.gembag)
    res$mean.recalls.gembag = cbind(rep(0,K),apply(res$recalls.gembag, 3, colMeans))
    res$mean.TPR.gembag =apply(res$TPR.gembag, 3, colMeans) 
    res$mean.FPR.gembag = apply(res$FPR.gembag, 3, colMeans) 
  }
  res$true.sparsity = spars.init
  return(res)
  
}

# Function for performing one iteration -----------------------------------------

jointGHS_simulation_one_iteration_ROC = function(n.vals,cov.matrices,prec.matrices,scale,include.jointGHS,include.SSJGL,include.JGL, include.gembag,
                                                 penalize.diagonal,seed,n.tau.jointGHS, n.lam.eff.ssjgl,n.lam1.jgl,n.thresh.gembag) {
  y = list()
  S_l = list()
  K=length(n.vals)
  p=ncol(prec.matrices[[1]])
  #tau_sq_vec = seq(1e-2,2,length.out=n.tau.jointGHS) 
  set.seed(seed)
  # Generate data. 
  for(k in 1:K){
    y[[k]] = mvtnorm::rmvnorm(n.vals[k], mean=rep(0,p), cov.matrices[[k]])
    S_l[[k]] = cov(y[[k]])
    if (scale) y[[k]] = scale(y[[k]])
  }
  # Perform joint methods
  if(include.jointGHS){
    #thresh.vals.jointGHS = c(seq(1e-7,0.01, length.out=n.tau.jointGHS/2), seq(0.01,0.4, length.out=n.tau.jointGHS/2))
    thresh.vals.jointGHS = c(seq(1e-9,0.01, length.out=n.tau.jointGHS/2), seq(0.01,0.4, length.out=n.tau.jointGHS/2))
    joint.theta = jointGHS(X=y, AIC_selection = T, AIC_eps = 0.1, epsilon=1e-3, verbose = F)$theta
    res.tmp = lapply(1:length(thresh.vals.jointGHS), FUN = function(i) lapply(1:K, FUN = function(k) cov2cor(joint.theta[[k]])>thresh.vals.jointGHS[i]))
  }
  else if(include.SSJGL){
    lam1 = 30
    lam2 = 1
    v1 = 1
    lam.eff = lam1 + seq(1,6,by=0.1)*5
    v0s = lam1/lam.eff
    ssjgl.tmp = SSJGL(Y=y,penalty='fused',lambda0=1, lambda1=lam1,lambda2=lam2, v1 = v1, v0s = v0s, tol.em=1e-4, a=1, b=p, doubly=TRUE, normalize=TRUE)
    ssjgl.tmp = ssjgl.tmp$thetalist # all estimates for different lam1.eff
  }
  else if(include.JGL){
    lambda1.vals = seq(0.4,1,length.out=n.lam1.jgl)
    jgl.tmp = lapply(1:length(lambda1.vals), FUN = function(i) JGL_select_AIC_lambda1fixed(Y=y,penalty='fused',lambda1= lambda1.vals[i],nlambda2=10,lambda2.min=0,lambda2.max=0.1,
                                                                                           penalize.diagonal=penalize.diagonal)$opt.fit)
  }
  else if(include.gembag){
    thresh.vals = seq(0.5,1,length.out=n.thresh.gembag) # For P
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
    # Threshold inclusion probability
    gembag.tmp = lapply(1:length(thresh.vals), FUN = function(i) lapply(1:K, FUN = function(k) gembag.tmp$P[[k]]>=thresh.vals[i]))
  }
  est=list()
  # Results from jointGHS
  if(include.jointGHS){
    est$precisions =  matrix(sapply(1:n.tau.jointGHS, FUN = function(s) sapply(1:K,FUN=function(k) precision(prec.matrices[[k]]!=0, abs(res.tmp[[s]][[k]])>1e-5))), 
                             ncol=n.tau.jointGHS, nrow = K, byrow=F)
    est$recalls = matrix(sapply(1:n.tau.jointGHS, FUN = function(s) sapply(1:K,FUN=function(k) recall(prec.matrices[[k]]!=0, abs(res.tmp[[s]][[k]])>1e-5 ))), 
                         ncol=n.tau.jointGHS, nrow = K, byrow=F) 
    est$TPR = matrix(sapply(1:n.tau.jointGHS, FUN = function(s) sapply(1:K,FUN=function(k) TPR(prec.matrices[[k]]!=0, abs(res.tmp[[s]][[k]])>1e-5 ))), 
                     ncol=n.tau.jointGHS, nrow = K, byrow=F) 
    est$FPR = matrix(sapply(1:n.tau.jointGHS, FUN = function(s) sapply(1:K,FUN=function(k) FPR(prec.matrices[[k]]!=0, abs(res.tmp[[s]][[k]])>1e-5 ))), 
                     ncol=n.tau.jointGHS, nrow = K, byrow=F) 
  }
  # Results from SSJGL
  else if(include.SSJGL){
    est$precisions.ssjgl =  matrix(sapply(1:n.lam.eff.ssjgl, FUN = function(s) sapply(1:K,FUN=function(k) precision(prec.matrices[[k]]!=0, ssjgl.tmp[[s]][[k]]!=0))), 
                                   ncol=n.lam.eff.ssjgl, nrow = K, byrow=F)
    est$recalls.ssjgl =  matrix(sapply(1:n.lam.eff.ssjgl, FUN = function(s) sapply(1:K,FUN=function(k) recall(prec.matrices[[k]]!=0, ssjgl.tmp[[s]][[k]]!=0))), 
                                ncol=n.lam.eff.ssjgl, nrow = K, byrow=F)
    est$TPR.ssjgl =  matrix(sapply(1:n.lam.eff.ssjgl, FUN = function(s) sapply(1:K,FUN=function(k) TPR(prec.matrices[[k]]!=0, ssjgl.tmp[[s]][[k]]!=0))), 
                            ncol=n.lam.eff.ssjgl, nrow = K, byrow=F)
    est$FPR.ssjgl =  matrix(sapply(1:n.lam.eff.ssjgl, FUN = function(s) sapply(1:K,FUN=function(k) FPR(prec.matrices[[k]]!=0, ssjgl.tmp[[s]][[k]]!=0))), 
                            ncol=n.lam.eff.ssjgl, nrow = K, byrow=F)
    
  }
  # Results from JGL
  else if(include.JGL){
    est$precisions.jgl = matrix(sapply(1:n.lam1.jgl, FUN = function(s) sapply(1:K,FUN=function(k) precision(prec.matrices[[k]]!=0, jgl.tmp[[s]][[k]]!=0))), 
                                ncol=n.lam1.jgl, nrow = K, byrow=F)
    est$recalls.jgl = matrix(sapply(1:n.lam1.jgl, FUN = function(s) sapply(1:K,FUN=function(k) recall(prec.matrices[[k]]!=0, jgl.tmp[[s]][[k]]!=0))), 
                             ncol=n.lam1.jgl, nrow = K, byrow=F)
    est$TPR.jgl = matrix(sapply(1:n.lam1.jgl, FUN = function(s) sapply(1:K,FUN=function(k) TPR(prec.matrices[[k]]!=0, jgl.tmp[[s]][[k]]!=0))), 
                         ncol=n.lam1.jgl, nrow = K, byrow=F)
    est$FPR.jgl = matrix(sapply(1:n.lam1.jgl, FUN = function(s) sapply(1:K,FUN=function(k) FPR(prec.matrices[[k]]!=0, jgl.tmp[[s]][[k]]!=0))), 
                         ncol=n.lam1.jgl, nrow = K, byrow=F)
  }
  # Results from GemBag
  else if(include.gembag){
    est$precisions.gembag = matrix(sapply(1:n.thresh.gembag, FUN = function(s) sapply(1:K,FUN=function(k) precision(prec.matrices[[k]]!=0, gembag.tmp[[s]][[k]]!=0))), 
                                   ncol=n.thresh.gembag, nrow = K, byrow=F)
    est$recalls.gembag = matrix(sapply(1:n.thresh.gembag, FUN = function(s) sapply(1:K,FUN=function(k) recall(prec.matrices[[k]]!=0, gembag.tmp[[s]][[k]]!=0))), 
                                ncol=n.thresh.gembag, nrow = K, byrow=F)
    est$TPR.gembag = matrix(sapply(1:n.thresh.gembag, FUN = function(s) sapply(1:K,FUN=function(k) TPR(prec.matrices[[k]]!=0, gembag.tmp[[s]][[k]]!=0))), 
                            ncol=n.thresh.gembag, nrow = K, byrow=F)
    est$FPR.gembag = matrix(sapply(1:n.thresh.gembag, FUN = function(s) sapply(1:K,FUN=function(k) FPR(prec.matrices[[k]]!=0, gembag.tmp[[s]][[k]]!=0))), 
                            ncol=n.thresh.gembag, nrow = K, byrow=F)
  }
  return(est)
}


JGL_select_AIC_lambda1fixed = function(Y,penalty='fused',lambda1,nlambda2,lambda2.min,lambda2.max, 
                                       penalize.diagonal){
  # JGL with lambda1 fixed, and lambda2 selecetd by the adapted AIC crierion (Danaher et al.)
  K=length(Y)
  p=ncol(Y[[1]])
  n.vals = unlist(lapply(Y,nrow))
  sample.cov = lapply(Y,cov)
  lambda2.vals = seq(lambda2.min,lambda2.max,length.out=nlambda2)
  mods.lam2=list()
  aic.lam2=rep(0,length(lambda2.vals))
  for (i in 1:length(lambda2.vals)){
    mod.temp = JGL(Y,penalty=penalty,lambda1=lambda1,lambda2 = lambda2.vals[i],penalize.diagonal = penalize.diagonal,
                   return.whole.theta = T)$theta
    mods.lam2[[i]] = mod.temp
    aic.lam2[i] = AIC_adapted(mod.temp,sample.cov=sample.cov,n.vals=n.vals)
  }
  opt.ind = which.min(aic.lam2) 
  res=list(opt.fit=mods.lam2[[opt.ind]],opt.lambda1 = lambda1,opt.lambda2 = lambda2.vals[opt.ind],
           opt.sparsities = unlist(lapply(mods.lam2[[opt.ind]],sparsity)))
  return(res)
}

TPR = function(g, g.hat){
  p <- nrow(g[, ])
  g <- g[, ]!=0
  g.hat <- g.hat[, ]!=0
  diag(g) <- rep(0, p)
  diag(g.hat[, ]) <- rep(0, p)
  tp = sum(g.hat[, ] == 1 & g[, ] == 1)/2
  pos = sum(g[, ] == 1)/2
  return(tp/pos)
}

FPR = function(g, g.hat){
  p <- nrow(g[, ])
  g <- g[, ]!=0
  g.hat <- g.hat[, ]!=0
  diag(g) <- rep(0, p)
  diag(g.hat[, ]) <- rep(0, p)
  fp = sum(g.hat[, ] == 1 & g[, ] == 0)/2
  neg = (sum(g[, ] == 0)-p)/2
  return(fp/neg)
}


