source('simulation_functions/perform_jointGHS_simulation.R')
source('simulation_functions/help_functions.R')
library(ggplot2)

# Test selection of tau with AIC, selecting the minimal value where the AIC has converged (for single-network version first)

# EXAMPLE 1 ------------------------------------------------------------------
# GENERATE GRAPH: n=150, p=100

n=150
p=100
set.seed(12345)
data.sf = huge::huge.generator(n=n, d=p,graph = 'scale-free') 
g.true.sf = data.sf$theta # True adjacency matrix
theta.true = data.sf$omega # The precision matrix
theta.true[which(theta.true<10e-5,arr.ind=T)]=0  
g.sf=graph.adjacency(data.sf$theta,mode="undirected",diag=F) # true igraph object
x.sf = data.sf$data # Observed attributes. nxp matrix.
x.sf.scaled= scale(x.sf) # Scale columns/variables.
s.sf.scaled = cov(x.sf.scaled) # Empirical covariance matrix
data.sf$sparsity
# 0.02


# Perform GHS for a range of tau_sq values

tau_sq_vals = seq(0.0001,1,length.out = 200)
sparsities.aic = rep(0,length(tau_sq_vals))
aic.vals = rep(0,length(tau_sq_vals))
diff.aic = 1
eps = 1e-1
i = 1
while(diff.aic > eps & i <= length(tau_sq_vals)){
  res = fastGHS::fastGHS(x.sf.scaled,tau_sq = tau_sq_vals[i],epsilon = 1e-3, fix_tau=TRUE)
  theta.est = res$theta
  theta.est[which(abs(theta.est)<1e-5, arr.ind=T)]=0
  prec.mat = cov2cor(theta.est)
  aic.vals[i] = JoStARS::gaussianAIC(s.sf.scaled,prec.mat, n)
  sparsities.aic[i] = tailoredGlasso::sparsity(theta.est!=0)
  if(i!=1){
    diff.aic = abs(aic.vals[i]-aic.vals[i-1])
  }
  i = i+1
}


# Optimal model 

ind.opt.aic = i-1
tau_sq_opt.aic =  tau_sq_vals[i-1]
tau_sq_opt.aic
# 0.1257156

mod.opt.aic = fastGHS::fastGHS(x.sf.scaled,tau_sq = tau_sq_opt.aic, epsilon=1e-3,fix_tau=TRUE)

theta.est.aic <- cov2cor(mod.opt.aic$theta)
theta.est.aic[which(abs(theta.est.aic) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.aic!=0)
# 0.006868687
tailoredGlasso::precision(as.matrix(theta.true!=0), theta.est.aic!=0)
# 0.7941176
tailoredGlasso::recall(as.matrix(theta.true!=0), theta.est.aic!=0)
# 0.2727273

# Works well!




