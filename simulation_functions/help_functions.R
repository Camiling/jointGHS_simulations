# Other functions -----------------------------------

# Find Matthews correlation coefficient for estimated graph
MCC = function(g,g.hat){
  p = nrow(g[,])
  diag(g) = rep(0,p) # Let diagonal elements be zero
  diag(g.hat) = rep(0,p) 
  tp = sum(g.hat ==1 & g ==1)/10 # True positives. Divide by 10 to avoid integer overflow. 
  fp = sum(g.hat ==1 & g ==0)/10 # False positives
  tn = (sum(g.hat == 0 & g == 0) - p)/10 # True negatives (do not include diagonal elements)
  fn = sum(g.hat == 0 & g == 1)/10 # False negatives
  return((tp*tn - fp*fn)/(sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))))
}



specificity = function(g,g.hat){
  # Specificity. How many of the edges that should have been negative are? TN/(TN+FP)
  conf.mat = confusion.matrix(g,g.hat)
  if(conf.mat[1,2]==0 & conf.mat[2,2]==0) return(1) # Avoid zero division
  else return(conf.mat[2,2]/(conf.mat[1,2]+conf.mat[2,2]))
}

FDR = function(g,g.hat){
  # False Discovery Rate
  return(1-precision(g,g.hat))
}

mutate.graph.alt = function(graph,fraction, generate.data=F, scale=T, larger.partialcor=F){
  # Mutate a given fraction of the edges of a graph. 
  # graph is the huge.generate() object to mutate, fraction is the fraction of edges to change. 
  # We basically 'swap pairs of nodes' by switching their cols and rows. 
  prec.mat = graph$omega
  prec.mat[which(abs(prec.mat)<10^(-11),arr.ind=T)]=0
  if(scale) cov.mat = cov2cor(graph$sigma) # added this
  else cov.mat = graph$sigma
  adj.mat = graph$theta
  data=graph$data
  p = ncol(graph$omega)
  
  adj.mat.upper = adj.mat
  adj.mat.upper[lower.tri(adj.mat.upper)]=0
  diag(adj.mat.upper) =0
  edges = which(as.matrix(adj.mat.upper)==1,arr.ind=T) # Edge pairs.
  n.mutations = floor(nrow(edges)*fraction)
  
  if(fraction==1){ # added this for jointGHS
    ans = list()
    if(larger.partialcor) graph.new = huge.generator(nrow(data),p,graph='scale-free',verbose = F,v=1,u=0.01)
    else graph.new = huge.generator(nrow(data),p,graph='scale-free',verbose = F,v=0.5,u=0.05)
    if(scale) ans$cov.mat = cov2cor(graph.new$sigma)
    else ans$cov.mat = graph.new$sigma
    ans$prec.mat = graph.new$omega
    ans$prec.mat[which(abs(ans$prec.mat)<10^(-11),arr.ind=T)]=0
    ans$adj.mat = graph.new$theta
    if(generate.data){
      ans$data = mvtnorm::rmvnorm(nrow(data), mean=rep(0,ncol(data)), ans$cov.mat)
    }
    return(ans)
  }
  
  if(n.mutations==0 | is.na(n.mutations)){
    ans = list()
    ans$cov.mat=cov.mat
    ans$prec.mat = prec.mat
    ans$adj.mat = adj.mat
    #ans$data = data # removed: should not reuse data
    if(generate.data){
      ans$data = mvtnorm::rmvnorm(nrow(data), mean=rep(0,ncol(data)), ans$cov.mat)
    }
    return(ans)
  }
  
  edges.to.change.ind = sample(1:nrow(edges),n.mutations) # We let the first index stay, then change the second one. 
  edges.to.change = edges[edges.to.change.ind,] # n.mutations x 2
  nodes.add = sample(1:p,n.mutations) # id of nodes to give the edges
  while(any(edges.to.change[,1] == nodes.add)){ # cannot give edge to itself
    edges.to.change.ind = sample(1:nrow(edges),n.mutations) # We let the first index stay, then change the second one. 
    edges.to.change = edges[edges.to.change.ind,] # n.mutations x 2
    nodes.add = sample(1:p,n.mutations) # id of nodes to give the edges
  }
  for(i in 1:n.mutations){
    tmp.prec=prec.mat
    tmp.adj = adj.mat
    tmp.dat = data
    tmp.cov.mat = cov.mat
    id.stay = edges.to.change[i,1]
    id.remove = edges.to.change[i,2]
    id.add=nodes.add[i]
    # Swap prec mat elements in rows. Then cols, the order does not matter!
    prec.mat[id.stay,id.add] = tmp.prec[id.stay,id.remove]
    prec.mat[id.stay,id.remove] = tmp.prec[id.stay,id.add]
    prec.mat[id.add,id.stay] = tmp.prec[id.remove,id.stay]
    prec.mat[id.remove,id.stay] = tmp.prec[id.add,id.stay]
    # swap adj mat rows
    adj.mat[id.stay,id.add] = tmp.adj[id.stay,id.remove]
    adj.mat[id.stay,id.remove] = tmp.adj[id.stay,id.add]
    adj.mat[id.add,id.stay] = tmp.adj[id.remove,id.stay]
    adj.mat[id.remove,id.stay] = tmp.adj[id.add,id.stay]
  }
  ans = list()
  if(scale) ans$cov.mat=cov2cor(solve(prec.mat))
  else ans$cov.mat= prec.mat
  ans$cov.mat[which(abs(ans$cov.mat)<1e-11,arr.ind = T)] = 0
  ans$prec.mat = prec.mat
  ans$adj.mat = adj.mat
  if(generate.data){
    # Generate new data
    ans$data = mvtnorm::rmvnorm(nrow(data), mean=rep(0,ncol(data)), ans$cov.mat)
  }
  return(ans)
}

mutate.graph= function(graph,fraction, generate.data=F, scale=T, larger.partialcor=F){
  # Mutate a given fraction of the edges of a graph. 
  # graph is the huge.generate() object to mutate, fraction is the fraction of edges to change. 
  # We basically 'swap pairs of nodes' by switching their cols and rows. 
  if(scale) prec.mat = cov2cor(graph$omega)
  else prec.mat = graph$omega
  prec.mat[which(abs(prec.mat)<10^(-4),arr.ind=T)]=0
  if(scale) cov.mat = cov2cor(graph$sigma) # added this
  else cov.mat = graph$sigma
  adj.mat = graph$theta
  data=graph$data
  p = ncol(graph$omega)
  
  adj.mat.upper = adj.mat
  adj.mat.upper[lower.tri(adj.mat.upper)]=0
  diag(adj.mat.upper) =0
  edges = which(as.matrix(adj.mat.upper)==1,arr.ind=T) # Edge pairs.
  n.mutations = floor(nrow(edges)*fraction)
  
  if(fraction==1){ # added this for jointGHS
    ans = list()
    if(larger.partialcor) graph.new = huge.generator(nrow(data),p,graph='scale-free',verbose = F,v=1,u=0.01)
    else graph.new = huge.generator(nrow(data),p,graph='scale-free',verbose = F,v=0.5,u=0.05)
    if(scale) ans$cov.mat = cov2cor(graph.new$sigma)
    else ans$cov.mat = graph.new$sigma
    if(scale) ans$prec.mat = cov2cor(graph.new$omega) 
    else ans$prec.mat = graph.new$omega
    ans$prec.mat[which(abs(ans$prec.mat)<10^(-4),arr.ind=T)]=0
    ans$adj.mat = graph.new$theta
    if(generate.data){
      ans$data = mvtnorm::rmvnorm(nrow(data), mean=rep(0,ncol(data)), ans$cov.mat)
    }
    return(ans)
  }
  
  if(n.mutations==0 | is.na(n.mutations)){
    ans = list()
    ans$cov.mat=cov.mat
    ans$prec.mat = prec.mat
    ans$adj.mat = adj.mat
    #ans$data = data # removed: should not reuse data
    if(generate.data){
      ans$data = mvtnorm::rmvnorm(nrow(data), mean=rep(0,ncol(data)), ans$cov.mat)
    }
    return(ans)
  }
  
  edges.to.change.ind = sample(1:nrow(edges),n.mutations) # We let the first index stay, then change the second one. 
  edges.to.change = edges[edges.to.change.ind,] # n.mutations x 2
  nodes.add = sample(1:p,n.mutations) # id of nodes to give the edges
  for(i in 1:n.mutations){
    tmp.prec=prec.mat
    tmp.adj = adj.mat
    tmp.dat = data
    tmp.cov.mat = cov.mat
    id.stay = edges.to.change[i,1]
    id.remove = edges.to.change[i,2]
    id.add=nodes.add[i]
    # Swap prec mat elements in rows. Then cols, the order does not matter!
    prec.mat[id.stay,id.add] = tmp.prec[id.stay,id.remove]
    prec.mat[id.stay,id.remove] = tmp.prec[id.stay,id.add]
    prec.mat[id.add,id.stay] = tmp.prec[id.remove,id.stay]
    prec.mat[id.remove,id.stay] = tmp.prec[id.add,id.stay]
    # swap adj mat rows
    adj.mat[id.stay,id.add] = tmp.adj[id.stay,id.remove]
    adj.mat[id.stay,id.remove] = tmp.adj[id.stay,id.add]
    adj.mat[id.add,id.stay] = tmp.adj[id.remove,id.stay]
    adj.mat[id.remove,id.stay] = tmp.adj[id.add,id.stay]
  }
  ans = list()
  ans$cov.mat=solve(prec.mat)
  ans$cov.mat[which(abs(ans$cov.mat)<1e-4,arr.ind = T)] = 0
  ans$prec.mat = prec.mat
  ans$adj.mat = adj.mat
  if(generate.data){
    # Generate new data
    ans$data = mvtnorm::rmvnorm(nrow(data), mean=rep(0,ncol(data)), ans$cov.mat)
  }
  return(ans)
}


JGL_select_AIC = function(Y,penalty='fused',nlambda1,lambda1.min,lambda1.max,nlambda2,lambda2.min,lambda2.max, 
                          lambda2.init,penalize.diagonal){
  # JGL with parameters selecetd by the adapted AIC crierion (Danaher et al.)
  K=length(Y)
  p=ncol(Y[[1]])
  n.vals = unlist(lapply(Y,nrow))
  sample.cov = lapply(Y,cov)
  lambda1.vals = seq(lambda1.min,lambda1.max,length.out=nlambda1)
  lambda2.vals = seq(lambda2.min,lambda2.max,length.out=nlambda2)
  mods.lam1=list()
  aic.lam1=rep(0,length(lambda1.vals)) 
  for (i in 1:length(lambda1.vals)){
    mod.temp = JGL(Y,penalty=penalty,lambda1=lambda1.vals[i],lambda2 = lambda2.init,return.whole.theta = T,penalize.diagonal=penalize.diagonal)$theta
    mods.lam1[[i]] = mod.temp
    aic.lam1[i] = AIC_adapted(mod.temp,sample.cov=sample.cov,n.vals=n.vals)
  }
  opt.ind.lam1 = which.min(aic.lam1)
  lambda1=lambda1.vals[opt.ind.lam1]
  
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

AIC_adapted = function(theta, sample.cov,n.vals) {
  if(length(theta) != length(sample.cov)) stop('number of precision matrices must be the same as number of covariance matrices')
  if (length(unique(c(unlist(lapply(theta, dim)),unlist(lapply(sample.cov, dim)))))!=1) stop("matrices must have the same dimension")
  if (any(unlist(lapply(theta, det)) <= 0 )) stop("precision matrix must be positive definite in all classes")
  if(any(unlist(lapply(sample.cov, isSymmetric)) != TRUE )) stop("sample covariance matrix must be symmetric")
  if (any(n.vals <= 0)) stop("number of observations n must be positive for all classes")
  p = dim(sample.cov[[1]])[2]
  K = length(sample.cov)
  aic.vals = rep(0,K)
  for(k in 1:length(theta))
  {
    theta.temp=theta[[k]]
    diag(theta.temp) = rep(0,p) # Do not count diagonal
    d=sum(theta.temp!=0)/2 # Number of edges
    aic.vals[k]  =  2*d - n.vals[k]*log(det(theta[[k]])) + n.vals[k]*sum(diag(sample.cov[[k]]%*%theta[[k]]))
  }
  return(sum(aic.vals))
}

matrix.distance <- function(mat1, mat2) {
  p <- nrow(mat1)
  if (mean(dim(mat1) == dim(mat2)) != 1) stop("matrices must have the same dimension")
  if (any(diag(mat1) == 0) | any(diag(mat2) == 0)) stop("diagonal elements of precision matrices cannot be zero")
  # Disregard almost-zero entries:
  mat1[which(abs(mat1) < 10^(-4), arr.ind = T)] <- 0
  mat2[which(abs(mat2) < 10^(-4), arr.ind = T)] <- 0
  m1 <- stats::cov2cor(as.matrix(Matrix::forceSymmetric(mat1)))
  m2 <- stats::cov2cor(as.matrix(Matrix::forceSymmetric(mat2)))
  observed.dist <- sum(abs(abs(m1) - abs(m2)))
  expected.dist <- (sum(abs(m1)) + sum(abs(m2)) - 2 * p)
  if (expected.dist == 0) { # No edges in either graph.
    mat.dist <- 0
  }
  else {
    mat.dist <- observed.dist / expected.dist
  }
  return(mat.dist)
}


print_results_jointGHS = function(obj.list,fracs.mutated, include.jointGHS=TRUE, include.GHS=TRUE, include.SSJGL=TRUE, include.JGL=TRUE, show.distance=F, show.interval=F, show.sd = F, 
                             show.specificity=F, collapse.values =F ){
  # obj is a list of objects returned by perform_jointGHS_simulation.
  # fracs.mutated is a vector of the mutated fraction in each simulation object
  # show.distance: should the matrix distance be printed?
  # show.interval: should intervals with the 2.5% and 97% quantiles be posted?
  # show.specificity: should the matrix distance be printed?
  # collapse.values: should the results for the different networks be collapsed?
  # Function for printing mean sparsities, precisions, recalls and matrix distances when several data sets were generated.
  # Note that we print the results for the different graphs on the same lines. 
  if(collapse.values){
    print_results_jointGHS_collapsed(obj.list,fracs.mutated,include.jointGHS=include.jointGHS, include.JGL=include.JGL, include.GHS=include.GHS, include.SSJGL=include.SSJGL,
                            show.distance=show.distance,show.specificity=show.specificity,show.interval=show.interval, show.sd = show.sd)
  }
  else if(show.sd==T){
    print_results_jointGHS_show_SD(obj.list,fracs.mutated=fracs.mutated,include.jointGHS=include.jointGHS, include.JGL=include.JGL,include.GHS=include.GHS, 
                                   include.SSJGL=include.SSJGL,show.distance=show.distance, show.specificity=show.specificity)
  }
  else{
    K = length(obj.list[[1]]$mean.opt.sparsities)
    if(show.interval){
      # Loop over each scenario
      for (i in 1:length(obj.list)){
        obj=obj.list[[i]]
        cat(fracs.mutated[i])
        if(include.GHS){
          cat(' & GHS ')
          for(k in 1:K){
            cat(' && ',round(obj$mean.opt.sparsities.ghs[k],3), '[',paste(round(quantile(obj$opt.sparsities.ghs[,k],probs=c(.025,.975)),3),collapse=','),'] &',
                round(obj$mean.precisions.ghs[k],2),'[',paste(round(quantile(obj$precisions.ghs[,k],probs=c(.025,.975)),3),collapse=','),'] &',
                round(obj$mean.recalls.ghs[k],2),'[',paste(round(quantile(obj$recalls.ghs[,k],probs=c(.025,.975)),3),collapse=','),']')
            if(show.specificity)cat('&',round(obj$mean.specificities.ghs[k],2), '[',paste(round(quantile(obj$specificities.ghs[,k],probs=c(.025,.975)),3),collapse=','),']')
            if(show.distance) cat(' & ',round(obj$mean.matrix.distances.ghs[k],3))
          }  
          cat(' \\\\ \n')
        }

        if(include.JGL){
          cat('  & JGL ')  
          for(k in 1:K){
            cat(' && ',round(obj$mean.opt.sparsities.jgl[k],3), '[',paste(round(quantile(obj$opt.sparsities.jgl[,k],probs=c(.025,.975)),3),collapse=','),'] &',
                round(obj$mean.precisions.jgl[k],2),'[',paste(round(quantile(obj$precisions.jgl[,k],probs=c(.025,.975)),3),collapse=','),'] &',
                round(obj$mean.recalls.jgl[k],2),'[',paste(round(quantile(obj$recalls.jgl[,k],probs=c(.025,.975)),3),collapse=','),']')
            if(show.specificity)cat('&',round(obj$mean.specificities.jgl[k],2),'[',paste(round(quantile(obj$specificities.jgl[,k],probs=c(.025,.975)),3),collapse=','),']') 
            if(show.distance) cat(' & ',round(obj$mean.matrix.distances.jgl[k],3))
          }  
          cat(' \\\\ \n')
        }
        if(include.SSJGL){
          cat('  & SSJGL ')  
          for(k in 1:K){
            cat(' && ',round(obj$mean.opt.sparsities.ssjgl[k],3), '[',paste(round(quantile(obj$opt.sparsities.ssjgl[,k],probs=c(.025,.975)),3),collapse=','),'] &',
                round(obj$mean.precisions.ssjgl[k],2),'[',paste(round(quantile(obj$precisions.ssjgl[,k],probs=c(.025,.975)),3),collapse=','),'] &',
                round(obj$mean.recalls.ssjgl[k],2),'[',paste(round(quantile(obj$recalls.ssjgl[,k],probs=c(.025,.975)),3),collapse=','),']')
            if(show.specificity)cat('&',round(obj$mean.specificities.ssjgl[k],2),'[',paste(round(quantile(obj$specificities.ssjgl[,k],probs=c(.025,.975)),3),collapse=','),']') 
            if(show.distance) cat(' & ',round(obj$mean.matrix.distances.ssjgl[k],3))
          }  
          cat(' \\\\ \n')
        }
        if(include.jointGHS){
          cat('  & jointGHS')
          for(k in 1:K){
            cat(' && ',round(obj$mean.opt.sparsities[k],3), '[',paste(round(quantile(obj$opt.sparsities[,k],probs=c(.025,.975)),3),collapse=','),'] &',
                round(obj$mean.precisions[k],2),'[',paste(round(quantile(obj$precisions[,k],probs=c(.025,.975)),3),collapse=','),'] &',
                round(obj$mean.recalls[k],2),'[',paste(round(quantile(obj$recalls[,k],probs=c(.025,.975)),3),collapse=','),']')
            if(show.specificity)cat('&',round(obj$mean.specificities[k],2), '[',paste(round(quantile(obj$specificities[,k],probs=c(.025,.975)),3),collapse=','),']') 
            if(show.distance) cat(' & ',round(obj$mean.matrix.distances[k],3))
          }  
          cat(' \\\\ \n \\hline \n')         
        }
      }
    }
    else{
      # Loop over each scenario
      for (i in 1:length(obj.list)){
        obj=obj.list[[i]]
        cat(fracs.mutated[i])
        if(include.GHS){
          cat(' & GHS ')
          for(k in 1:K){
            cat(' && ',round(obj$mean.opt.sparsities.ghs[k],3), '& ',
                round(obj$mean.precisions.ghs[k],2),' & ',
                round(obj$mean.recalls.ghs[k],2))
            if(show.specificity)cat('&',round(obj$mean.specificities.ghs[k],2))
            if(show.distance) cat(' & ',round(obj$mean.matrix.distances.ghs[k],3))
          }  
          cat(' \\\\ \n')
        }
        if(include.JGL){
          cat('  & JGL ')  
          for(k in 1:K){
            cat(' && ',round(obj$mean.opt.sparsities.jgl[k],3), '& ',
                round(obj$mean.precisions.jgl[k],2),'& ',
                round(obj$mean.recalls.jgl[k],2))
            if(show.specificity)cat('&',round(obj$mean.specificities.jgl[k],2)) 
            if(show.distance) cat(' & ',round(obj$mean.matrix.distances.jgl[k],3))
          }  
          cat(' \\\\ \n')
        }
        if(include.SSJGL){
          cat('  & SSJGL ')  
          for(k in 1:K){
            cat(' && ',round(obj$mean.opt.sparsities.ssjgl[k],3), '& ',
                round(obj$mean.precisions.ssjgl[k],2),'& ',
                round(obj$mean.recalls.ssjgl[k],2))
            if(show.specificity)cat('&',round(obj$mean.specificities.ssjgl[k],2)) 
            if(show.distance) cat(' & ',round(obj$mean.matrix.distances.ssjgl[k],3))
          }  
          cat(' \\\\ \n')
        }
        if(include.jointGHS){
          cat('  & jointGHS')
          for(k in 1:K){
            cat(' && ',round(obj$mean.opt.sparsities[k],3), '& ',
                round(obj$mean.precisions[k],2),' & ',
                round(obj$mean.recalls[k],2))
            if(show.specificity)cat('&',round(obj$mean.specificities[k],2)) 
            if(show.distance) cat(' & ',round(obj$mean.matrix.distances[k],3))
          }  
        }
        cat(' \\\\ \n \\hline \n')  
      }
    }
  }
}

print_results_jointGHS_show_SD = function(obj.list,fracs.mutated, include.jointGHS=TRUE, include.GHS=TRUE, include.SSJGL=TRUE, include.JGL=TRUE,show.distance=F,show.specificity=F){
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
    if(include.GHS){
      cat(' & GHS ')
      for(k in 1:K){
        cat(' && ',round(obj$mean.opt.sparsities.ghs[k],3), '(',round(sd(obj$opt.sparsities.ghs[,k]),3),')',' & ',
            round(obj$mean.precisions.ghs[k],2),'(',round(sd(obj$precisions.ghs[,k]),2),')',' & ',
            round(obj$mean.recalls.ghs[k],2), '(',round(sd(obj$recalls.ghs[,k]),2),')')
        if(show.specificity)cat('&',round(obj$mean.specificities.ghs[k],2), '(',round(sd(obj$specificities.ghs[,k]),2),')')
        if(show.distance) cat(' & ',round(obj$mean.matrix.distances.ghs[k],3))
      }
      cat(' \\\\ \n')
    }
    if(include.JGL){
      cat('  & JGL ')  
      for(k in 1:K){
        cat(' && ',round(obj$mean.opt.sparsities.jgl[k],3), '(',round(sd(obj$opt.sparsities.jgl[,k]),3),')',' & ',
            round(obj$mean.precisions.jgl[k],2),'(',round(sd(obj$precisions.jgl[,k]),2),')',' & ',
            round(obj$mean.recalls.jgl[k],2),'(',round(sd(obj$recalls.jgl[,k]),2),')')
        if(show.specificity)cat('&',round(obj$mean.specificities.jgl[k],2), '(',round(sd(obj$specificities.jgl[,k]),2),')') 
        if(show.distance) cat(' & ',round(obj$mean.matrix.distances.jgl[k],3))
      }  
      cat(' \\\\ \n')
    }
    if(include.SSJGL){
      cat('  & SSJGL ')  
      for(k in 1:K){
        cat(' && ',round(obj$mean.opt.sparsities.ssjgl[k],3), '(',round(sd(obj$opt.sparsities.ssjgl[,k]),3),')',' & ',
            round(obj$mean.precisions.ssjgl[k],2),'(',round(sd(obj$precisions.ssjgl[,k]),2),')',' & ',
            round(obj$mean.recalls.ssjgl[k],2),'(',round(sd(obj$recalls.ssjgl[,k]),2),')')
        if(show.specificity)cat('&',round(obj$mean.specificities.ssjgl[k],2), '(',round(sd(obj$specificities.ssjgl[,k]),2),')') 
        if(show.distance) cat(' & ',round(obj$mean.matrix.distances.ssjgl[k],3))
      }  
      cat(' \\\\ \n')
    }
    if(include.jointGHS){
      cat('  & jointGHS')
      for(k in 1:K){
        cat(' && ',round(obj$mean.opt.sparsities[k],3), '(',round(sd(obj$opt.sparsities[,k]),3),')',' & ',
            round(obj$mean.precisions[k],2),'(',round(sd(obj$precisions[,k]),2),')',' & ',
            round(obj$mean.recalls[k],2),'(',round(sd(obj$recalls[,k]),2),')')
        if(show.specificity)cat('&',round(obj$mean.specificities[k],2), '(',round(sd(obj$specificities[,k]),2),')') 
        if(show.distance) cat(' & ',round(obj$mean.matrix.distances[k],3))
      }  
    }
    cat(' \\\\ \n \\hline \n')
  }
}

print_results_jointGHS_collapsed = function(obj.list,fracs.mutated,include.jointGHS=TRUE, include.GHS=TRUE, include.SSJGL=TRUE, include.JGL=TRUE,show.distance=F,show.specificity=F,show.interval=F, show.sd =T){
  K = length(obj.list[[1]]$mean.opt.sparsities)
  
  # We find the average means and average sd. 
  # Intervals are not implemented, as they cannot be averaged
  # Loop over each scenario
  for (i in 1:length(obj.list)){
    obj=obj.list[[i]]
    cat(fracs.mutated[i])
    # Collapse values
    
    if(include.GHS){
      cat(' & GHS ')
      if(show.sd){
        sep1 = c('(',round(mean(apply(obj$opt.sparsities.ghs, 2, sd)),3),') & ')
        sep2 = c('(',round(mean(apply(obj$precisions.ghs, 2, sd)),2),') & ')
        sep3 = c('(',round(mean(apply(obj$recalls.ghs, 2, sd)),2),')')
        if (show.specificity) {
          sep4 = c('(',round(mean(apply(obj$specificities.ghs, 2, sd)),2),')')
        }
      }
      else{
        sep1 = sep2 = ' & '
        sep3  = sep4 = ""
      }
      cat(' && ',round(mean(colMeans(obj$opt.sparsities.ghs)),3), sep1,
          round(mean(colMeans(obj$precisions.ghs)),2), sep2,
          round(mean(colMeans(obj$recalls.ghs)),2), sep3)
      if(show.specificity)cat('&',round(mean(colMeans(obj$specificities.ghs)),2), sep4)
      if(show.distance) cat(' & ',round(mean(colMeans(obj$matrix.distances.ghs)),3))
      cat(' \\\\ \n')   
    }
    if(include.JGL){
      cat('  & JGL ')  
      if(show.sd){
        sep1 = c('(',round(mean(apply(obj$opt.sparsities.jgl, 2, sd)),3),') & ')
        sep2 = c('(',round(mean(apply(obj$precisions.jgl, 2, sd)),2),') & ')
        sep3 = c('(',round(mean(apply(obj$recalls.jgl, 2, sd)),2),')')
        if (show.specificity) {
          sep4 = c('(',round(mean(apply(obj$specificities.jgl, 2, sd)),2),')')
        }
      }
      else{
        sep1 = sep2 = ' & '
        sep3  = sep4 = ""
      }
      cat(' && ',round(mean(colMeans(obj$opt.sparsities.jgl)),3), sep1,
          round(mean(colMeans(obj$precisions.jgl)),2), sep2,
          round(mean(colMeans(obj$recalls.jgl)),2), sep3)
      if(show.specificity)cat('&',round(mean(colMeans(obj$specificities.jgl)),2), sep4) 
      if(show.distance) cat(' & ',round(mean(colMeans(obj$matrix.distances.jgl)),3))
      cat(' \\\\ \n')
    }
    if(include.SSJGL){
      cat('  & SSJGL ')  
      if(show.sd){
        sep1 = c('(',round(mean(apply(obj$opt.sparsities.ssjgl, 2, sd)),3),') & ')
        sep2 = c('(',round(mean(apply(obj$precisions.ssjgl, 2, sd)),2),') & ')
        sep3 = c('(',round(mean(apply(obj$recalls.ssjgl, 2, sd)),2),')')
        if (show.specificity) {
          sep4 = c('(',round(mean(apply(obj$specificities.ssjgl, 2, sd)),2),')')
        }
      }
      else{
        sep1 = sep2 = ' & '
        sep3  = sep4 = ""
      }
      cat(' && ',round(mean(colMeans(obj$opt.sparsities.ssjgl)),3), sep1,
          round(mean(colMeans(obj$precisions.ssjgl)),2), sep2,
          round(mean(colMeans(obj$recalls.ssjgl)),2), sep3)
      if(show.specificity)cat('&',round(mean(colMeans(obj$specificities.ssjgl)),2), sep4) 
      if(show.distance) cat(' & ',round(mean(colMeans(obj$matrix.distances.ssjgl)),3))
      cat(' \\\\ \n')
    }
    if(include.jointGHS){
      cat('  & jointGHS')
      if(show.sd){
        sep1 = c('(',round(mean(apply(obj$opt.sparsities, 2, sd)),3),') & ')
        sep2 = c('(',round(mean(apply(obj$precisions, 2, sd)),2),') & ')
        sep3 = c('(',round(mean(apply(obj$recalls, 2, sd)),2),')')
        if (show.specificity) {
          sep4 = c('(',round(mean(apply(obj$specificities, 2, sd)),2),')')
        }
      }
      else{
        sep1 = sep2 = ' & '
        sep3  = sep4 = ""
      }
      cat(' && ',round(mean(colMeans(obj$opt.sparsities)),3),sep1,
          round(mean(colMeans(obj$precisions)),2),sep2,
          round(mean(colMeans(obj$recalls)),2), sep3)
      if(show.specificity)cat('&',round(mean(colMeans(obj$specificities)),2), sep4) 
      if(show.distance) cat(' & ',round(mean(colMeans(obj$matrix.distances)),3))
    }
    cat(' \\\\ \n \\hline \n')
  }
}





print_results_fastGHS = function(obj.list,include.glasso=TRUE, include.GHS=TRUE, include.fastGHS=TRUE, show.distance=F,show.interval=T, show.sd = F, 
                                  show.specificity=F){
  # obj is a list of objects returned by perform_jointGHS_simulation.
  # show.distance: should the matrix distance be printed?
  # show.interval: should intervals with the 2.5% and 97% quantiles be posted?
  # show.specificity: should the matrix distance be printed?
  # Function for printing mean sparsities, precisions, recalls and matrix distances when several data sets were generated.
  # Note that we print the results for the different graphs on the same lines. 
  if(show.sd==T){
    print_results_fastGHS_show_SD(obj.list,include.glasso=include.glasso, include.GHS=include.GHS, include.fastGHS=include.fastGHS, show.distance, show.specificity)
  }
  else{
    if(show.interval){
      # Loop over each scenario
      for (i in 1:length(obj.list)){
        obj=obj.list[[i]]
        if(include.glasso){
          cat(' & Glasso ')
          cat(' && ',round(obj$mean.opt.sparsities.glasso,3),'[',paste(round(quantile(obj$opt.sparsities.glasso,probs=c(.025,.975)),3),collapse=','),'] &',
              round(obj$mean.precisions.glasso,2),'[',paste(round(quantile(obj$precisions.glasso,probs=c(.025,.975)),3),collapse=','),'] &',
              round(obj$mean.recalls.glasso,2),'[',paste(round(quantile(obj$recalls.glasso,probs=c(.025,.975)),3),collapse=','),']')
          if(show.specificity)cat('&',round(obj$mean.specificities.glasso,2),'[',paste(round(quantile(obj$specificities.glasso,probs=c(.025,.975)),3),collapse=','),']')
          if(show.distance) cat(' & ',round(obj$mean.matrix.distances.glasso,3))
          cat(' \\\\ \n')         
        }
        if(include.GHS){
          cat(' & GHS ')
          cat(' && ',round(obj$mean.opt.sparsities.ghs,3), '[',paste(round(quantile(obj$opt.sparsities.ghs,probs=c(.025,.975)),3),collapse=','),'] &',
              round(obj$mean.precisions.ghs,2),'[',paste(round(quantile(obj$precisions.ghs,probs=c(.025,.975)),3),collapse=','),'] &',
              round(obj$mean.recalls.ghs,2),'[',paste(round(quantile(obj$recalls.ghs,probs=c(.025,.975)),3),collapse=','),']')
          if(show.specificity)cat('&',round(obj$mean.specificities.ghs,2), '[',paste(round(quantile(obj$specificities.ghs,probs=c(.025,.975)),3),collapse=','),']')
          if(show.distance) cat(' & ',round(obj$mean.matrix.distances.ghs,3))
          cat(' \\\\ \n')          
        }
        if(include.fastGHS){
          cat(' & fastGHS ')
          cat(' && ',round(obj$mean.opt.sparsities.fastghs,3), '[',paste(round(quantile(obj$opt.sparsities.fastghs,probs=c(.025,.975)),3),collapse=','),'] &',
              round(obj$mean.precisions.fastghs,2),'[',paste(round(quantile(obj$precisions.fastghs,probs=c(.025,.975)),3),collapse=','),'] &',
              round(obj$mean.recalls.fastghs,2),'[',paste(round(quantile(obj$recalls.fastghs,probs=c(.025,.975)),3),collapse=','),']')
          if(show.specificity)cat('&',round(obj$mean.specificities.fastghs,2), '[',paste(round(quantile(obj$specificities.fastghs,probs=c(.025,.975)),3),collapse=','),']')
          if(show.distance) cat(' & ',round(obj$mean.matrix.distances.fastghs,3))
          cat(' \\\\ \n')          
        }
        cat(' \\\\ \n \\hline \n')
      }
    }
    else{
      # Loop over each scenario
      for (i in 1:length(obj.list)){
        obj=obj.list[[i]]
        cat(' & Glasso ')
        if(include.glasso){
          cat(' && ',round(obj$mean.opt.sparsities.glasso,3), ' & ',
              round(obj$mean.precisions.glasso,2),' & ',
              round(obj$mean.recalls.glasso,2))
          if(show.specificity)cat('&',round(obj$mean.specificities.glasso,2))
          if(show.distance) cat(' & ',round(obj$mean.matrix.distances.glasso,3))
          cat(' \\\\ \n')
        }    
        if(include.GHS){
          cat(' & GHS ')
          cat(' && ',round(obj$mean.opt.sparsities.ghs,3), '& ',
              round(obj$mean.precisions.ghs,2),' & ',
              round(obj$mean.recalls.ghs,2))
          if(show.specificity)cat('&',round(obj$mean.specificities.ghs,2))
          if(show.distance) cat(' & ',round(obj$mean.matrix.distances.ghs,3))
          cat(' \\\\ \n')          
        }
        if(include.fastGHS){
          cat('  & fastGHS')
          cat(' && ',round(obj$mean.opt.sparsities.fastghs,3), '& ',
              round(obj$mean.precisions.fastghs,2),' & ',
              round(obj$mean.recalls.fastghs,2))
          if(show.specificity)cat('&',round(obj$mean.specificities.fastghs,2)) 
          if(show.distance) cat(' & ',round(obj$mean.matrix.distances.fastghs,3))
          cat(' \\\\ \n \\hline \n')
        }

      }
    }
  }
}

print_results_fastGHS_show_SD = function(obj.list,include.glasso=TRUE, include.GHS=TRUE, include.fastGHS=TRUE, show.distance=F,show.specificity=F){
  # obj is a list of objects returned by perform_jointGHS_simulation.
  # fracs.mutated is a vector of the mutated fraction in each simulation object
  # show.distance: should the matrix distance be printed?
  # show.specificity: should the matrix distance be printed?
  # Function for printing mean sparsities, precisions, recalls and matrix distances when several data sets were generated.
  # Note that we print the results for the different graphs on the same lines. 
  # Loop over each scenario
  for (i in 1:length(obj.list)){
    obj=obj.list[[i]]
    if(include.glasso){
      cat(' & Glasso ')
        cat(' && ',round(obj$mean.opt.sparsities.glasso,3), '(',round(sd(obj$opt.sparsities.glasso),3),')',' & ',
            round(obj$mean.precisions.glasso,2),'(',round(sd(obj$precisions.glasso),2),')',' & ',
            round(obj$mean.recalls.glasso,2), '(',round(sd(obj$recalls.glasso),2),')')
        if(show.specificity)cat('&',round(obj$mean.specificities.glasso,2), '(',round(sd(obj$specificities.glasso),2),')')
        if(show.distance) cat(' & ',round(obj$mean.matrix.distances.glasso,3))
        cat(' \\\\ \n')
    }
    if(include.GHS){
      cat(' & GHS ')
      cat(' && ',round(obj$mean.opt.sparsities.ghs,3), '(',round(sd(obj$opt.sparsities.ghs),3),')',' & ',
          round(obj$mean.precisions.ghs,2),'(',round(sd(obj$precisions.ghs),2),')',' & ',
          round(obj$mean.recalls.ghs,2), '(',round(sd(obj$recalls.ghs),2),')')
      if(show.specificity)cat('&',round(obj$mean.specificities.ghs,2), '(',round(sd(obj$specificities.ghs),2),')')
      if(show.distance) cat(' & ',round(obj$mean.matrix.distances.ghs,3))
      cat(' \\\\ \n')
    }
    if(include.fastGHS){
      cat('  & fastGHS')
      cat(' && ',round(obj$mean.opt.sparsities.fastghs,3), '(',round(sd(obj$opt.sparsities.fastghs),3),')',' & ',
          round(obj$mean.precisions.fastghs,2),'(',round(sd(obj$precisions.fastghs),2),')',' & ',
          round(obj$mean.recalls.fastghs,2),'(',round(sd(obj$recalls.fastghs),2),')')
      if(show.specificity)cat('&',round(obj$mean.specificities.fastghs,2), '(',round(sd(obj$specificities.fastghs),2),')') 
      if(show.distance) cat(' & ',round(obj$mean.matrix.distances.fastghs,3))
    cat(' \\\\ \n \\hline \n')
    }
  }
}



gaussianAIC <- function(sample.cov, theta, n) {
  p <- nrow(theta)
  theta2 <- theta
  diag(theta2) <- rep(0, p)
  d <- sum(theta2 != 0) / 2
  return(-2 * tailoredGlasso::gaussianloglik(sample.cov, theta, n) + 2 * d)
}

