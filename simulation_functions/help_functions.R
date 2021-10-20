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


mutate.graph= function(graph,fraction, generate.data=F){
  # Mutate a given fraction of the edges of a graph. 
  # graph is the huge.generate() object to mutate, fraction is the fraction of edges to change. 
  # We basically 'swap pairs of nodes' by switching their cols and rows. 
  prec.mat = cov2cor(graph$omega) # added this for jointGHS. Include scale argument?
  prec.mat[which(abs(prec.mat)<10^(-4),arr.ind=T)]=0
  cov.mat = cov2cor(graph$sigma) # added this
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
    graph.new = huge.generator(nrow(data),p,graph='scale-free',verbose = F,v=0.5,u=0.05)
    ans$cov.mat = cov2cor(graph.new$sigma)
    ans$prec.mat = cov2cor(graph.new$omega) # added this for jointGHS. Include scale argument?
    ans$prec.mat[which(abs(ans$prec.mat)<10^(-4),arr.ind=T)]=0
    ans$adj.mat = graph.new$theta
    ans$data = mvtnorm::rmvnorm(nrow(data), mean=rep(0,ncol(data)), ans$cov.mat)
    return(ans)
  }
  
  if(n.mutations==0 | is.na(n.mutations)){
    ans = list()
    ans$cov.mat=cov.mat
    ans$prec.mat = prec.mat
    ans$adj.mat = adj.mat
    #ans$data = data # removed: should not reuse data
    ans$data = mvtnorm::rmvnorm(nrow(data), mean=rep(0,ncol(data)), ans$cov.mat)
    return(ans)
  }
  
  edges.to.change.ind = sample(1:nrow(edges),n.mutations) # We let the first index stay, then change the second one. 
  edges.to.change = edges[edges.to.change.ind,] # n.mutations x 2
  nodes.add = sample(1:p,n.mutations) # id of nodes to give the edges
  nodes.remove = edges[sample(1:nrow(edges),n.mutations),1] # The nodes to 'leave out'
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




