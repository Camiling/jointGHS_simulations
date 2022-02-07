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


