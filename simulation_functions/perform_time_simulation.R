#' Perform fastGHS simulation to assess time use
#' 
#' This function performs simulations for fastGHS, assessing time use 
#' 
#' @return simulation results, including sparsity, precision, recall and specificity
perform_time_simulation = function(p,n,nCores=5, AIC_selection=F, include.GHS=T){
  n.p = length(p)
  # Create data sets
  data.sets = lapply(p,FUN=function(pp) scale(huge.generator(n,pp,graph='scale-free',verbose = F)$data))
  # Perform fastGHS
  registerDoParallel(nCores)
  res.list = foreach (i=1:n.p) %dopar% {
    system.time(fastGHS::fastGHS(X=data.sets[[i]], AIC_selection = AIC_selection, epsilon = 1e-3, fix_tau = T,tau_sq=1, verbose=FALSE))['elapsed'];
  }
  registerDoSEQ()
  
  if(!include.GHS){
    return(unlist(res.list))
  }
  else{
    # Perform GHS when feasible
    if(any(p>=100)){
      p.ghs = p[1:(min(which(p>=100))-1)] # Only smaller than 100
    }
    else{
     p.ghs = p
    }
    n.p.ghs = length(p.ghs)
    res.ghs=rep(10000,n.p.ghs)
    for(i in 1:n.p.ghs){
      res.ghs[i] = system.time(GHS(t(data.sets[[i]])%*%(data.sets[[i]]),n,burnin=100,nmc=1000))['elapsed']
    }
  
    times.fastGHS = unlist(res.list)
    times.GHS = res.ghs
    return(list(times.fastGHS, times.GHS))
  }
}

#' Perform fastGHS simulation to assess time use
#' 
#' This function performs simulations for fastGHS, assessing time use 
#' 
#' @return simulation results, including sparsity, precision, recall and specificity
perform_time_simulation_withglasso = function(p,n,nCores=5, AIC_selection=F, include.GHS=T){
  n.p = length(p)
  # Create data sets
  data.sets = lapply(p,FUN=function(pp) scale(huge.generator(n,pp,graph='scale-free',verbose = F)$data))
  # Perform fastGHS
  registerDoParallel(nCores)
  res.list = foreach (i=1:n.p) %dopar% {
    system.time(fastGHS::fastGHS(X=data.sets[[i]], AIC_selection = AIC_selection, epsilon = 1e-3, fix_tau = T,tau_sq=1, verbose=FALSE))['elapsed'];
  }
  registerDoSEQ()
  
  res.list.glasso = foreach (i=1:n.p) %dopar% {
    system.time(huge.select(huge(data.sets[[i]],method='glasso',verbose = F),criterion='stars', stars.thresh = 0.03, verbose = F))['elapsed'];
  }
  
  if(!include.GHS){
    return(list(unlist(res.list), unlist(res.list.glasso)))
  }
  else{
    # Perform GHS 
    p.ghs = p
    n.p.ghs = length(p.ghs)
    res.ghs=rep(10000,n.p.ghs)
    for(i in 1:n.p.ghs){
      res.ghs[i] = system.time(GHS(t(data.sets[[i]])%*%(data.sets[[i]]),n,burnin=100,nmc=1000))['elapsed']
    }
    times.fastGHS = unlist(res.list)
    times.glasso = unlist(res.list.glasso)
    times.GHS = res.ghs
    return(list(times.fastGHS, times.GHS, times.glasso))
  }
}

GemBag_oneit = function(X_l,n.vals,p){
  # Perform joint method
  # Set of hyperparameters
  S_l = lapply(X_l, cov)
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
  return(gembag.tmp)
}

#' Perform joint network simulation to assess time use
#' 
#' This function performs simulations, assessing time use for all joint methods
#' 
#' @return simulation results, including sparsity, precision, recall and specificity
perform_time_simulation_joint = function(p,K,n.vals,nCores=5){
  n.p = length(p)
  # Create data sets
  data.sets=list()
  for(i in 1:length(p)){
    data.tmp = list()
    huge.init = huge.generator(n.vals[1],p[i],graph='scale-free',verbose = F,v=1,u=0.01)
    data.tmp[[1]] = scale(huge.init$data)
    for(k in 2:K){
      valid=F
      while(!valid){ # Ensure valid precision matrices
        huge.tmp = mutate.graph(huge.init,fraction=0.5,scale=T, generate.data = F)
        cov.mat = huge.tmp$cov.mat
        # Avoid rounding errors leading to matrices not being symmetric
        if(!matrixcalc::is.symmetric.matrix(cov.mat)){
          cov.mat = round(cov.mat,8)
        }
        valid = matrixcalc::is.positive.definite(cov.mat)
      }
      data.tmp[[k]] = scale(mvtnorm::rmvnorm(n.vals[k], mean=rep(0,p[i]), cov.mat))
    }
    data.sets[[i]] = data.tmp
  }
  # Perform jointGHS
  registerDoParallel(nCores)
  res.list = foreach (i=1:n.p) %dopar% {
    system.time(jointGHS::jointGHS(X=data.sets[[i]], AIC_selection = F,fix_tau = T, tau_sq= 1, epsilon = 1e-3 , verbose=FALSE))['elapsed'];
  }
  registerDoSEQ()
  times.jointGHS = unlist(res.list)
  cat('jointGHS done \n')
  # times.jointGHS = rep(0,n.p)
  # for(i in 1:n.p){
  #   times.jointGHS[i] = system.time(jointGHS::jointGHS(X=data.sets[[i]], AIC_selection = T, AIC_eps =5, tau_sq_min = 0.5,
  #                                                tau_sq_stepsize = 0.2, epsilon = 1e-3 ,nCores = nCores, verbose=FALSE))['elapsed'];
  # }
  # Perform SSJGL
  registerDoParallel(nCores)
  res.list = foreach (i=1:n.p) %dopar% {
    system.time(SSJGL(Y=data.sets[[i]],penalty='fused',lambda0=1, lambda1=1,lambda2=1, v1 = 1, v0s = 1/(1 + c(1:10) * 5), tol.em=1e-4, a=1, b=p[i],
                      doubly=TRUE, normalize=TRUE))['elapsed'];
  }
  registerDoSEQ()
  times.SSJGL = unlist(res.list)
  cat('SSJGL done \n')
  # Perform JGL
  registerDoParallel(nCores)
  res.list = foreach (i=1:n.p) %dopar% {
    system.time(JGL_select_AIC(Y=data.sets[[i]],penalty='fused',nlambda1=20,lambda1.min=0.001,lambda1.max=0.3,nlambda2=20,lambda2.min=0,lambda2.max=0.1,lambda2.init=0.01,
                               penalize.diagonal=FALSE))['elapsed'];
  }
  registerDoSEQ()
  times.JGL = unlist(res.list)
  cat('JGL done \n')
  # Perform GemBag
  registerDoParallel(nCores)
  res.list = foreach (i=1:n.p) %dopar% {
    system.time(GemBag_oneit(data.sets[[i]],n.vals,p[i]))['elapsed'];
  }
  registerDoSEQ()
  times.GemBag = unlist(res.list)
  
  return(list(times.jointGHS=times.jointGHS, times.SSJGL=times.SSJGL, times.JGL = times.JGL, times.GemBag = times.GemBag))
}

#' Perform joint network simulation to assess time use
#' 
#' This function performs simulations, assessing time use for all joint methods
#' 
#' @return simulation results, including sparsity, precision, recall and specificity
perform_time_simulation_joint_onlyjointGHS = function(p,K,n.vals,nCores=5){
  n.p = length(p)
  # Create data sets
  data.sets=list()
  for(i in 1:length(p)){
    data.tmp = list()
    huge.init = huge.generator(n.vals[1],p[i],graph='scale-free',verbose = F,v=1,u=0.01)
    data.tmp[[1]] = scale(huge.init$data)
    for(k in 2:K){
      valid=F
      while(!valid){ # Ensure valid precision matrices
        huge.tmp = mutate.graph(huge.init,fraction=0.5,scale=T, generate.data = F)
        cov.mat = huge.tmp$cov.mat
        # Avoid rounding errors leading to matrices not being symmetric
        if(!matrixcalc::is.symmetric.matrix(cov.mat)){
          cov.mat = round(cov.mat,8)
        }
        valid = matrixcalc::is.positive.definite(cov.mat)
      }
      data.tmp[[k]] = scale(mvtnorm::rmvnorm(n.vals[k], mean=rep(0,p[i]), cov.mat))
    }
    data.sets[[i]] = data.tmp
  }
  # Perform jointGHS
  registerDoParallel(nCores)
  res.list = foreach (i=1:n.p) %dopar% {
    system.time(jointGHS::jointGHS(X=data.sets[[i]], AIC_selection = F,fix_tau = T, tau_sq= 1, epsilon = 1e-3 , verbose=FALSE))['elapsed'];
  }
  registerDoSEQ()
  times.jointGHS = unlist(res.list)
  return(times.jointGHS)
}


