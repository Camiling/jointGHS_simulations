source('simulation_functions/perform_jointGHS_simulation.R')
source('simulation_functions/help_functions.R')
library(ggplot2)

# Test selection of tau by eBIC and AIC (for single-network version first)

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

tau_sq_vals = seq(0.0001,2,length.out = 100)
ebic.vals = rep(0,length(tau_sq_vals))
sparsities.ebic = rep(0,length(tau_sq_vals))
aic.vals = rep(0,length(tau_sq_vals))
logliks.ebic = rep(0, length(tau_sq_vals))
for(i in 1:length(tau_sq_vals)){
  res <- fastGHS::fastGHS(x.sf.scaled,tau_sq = tau_sq_vals[i],epsilon = 1e-4, fix_tau=TRUE)
  theta.est = res$theta
  theta.est[which(abs(theta.est)<1e-5, arr.ind=T)]=0
  prec.mat = cov2cor(theta.est)
  ebic.vals[i] = JoStARS::eBIC(s.sf.scaled,prec.mat, n, gamma=0)
  aic.vals[i] = JoStARS::gaussianAIC(s.sf.scaled,prec.mat, n)
  sparsities.ebic[i] = tailoredGlasso::sparsity(theta.est!=0)
  logliks.ebic[i] = tailoredGlasso::gaussianloglik(s.sf.scaled,prec.mat, n)
}

# Optimal model found by eBIC
ind.opt = which.min(ebic.vals)
tau_sq_opt = tau_sq_vals[ind.opt]
mod.opt = fastGHS::fastGHS(x.sf.scaled,tau_sq = tau_sq_opt, epsilon=1e-4,fix_tau=TRUE)

theta.est <- cov2cor(mod.opt$theta)
theta.est[which(abs(theta.est) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est!=0)
# 0
tailoredGlasso::precision(as.matrix(theta.true!=0), theta.est!=0)
# 0
tailoredGlasso::recall(as.matrix(theta.true!=0), theta.est!=0)
# 0

tau_sq_opt

# Selects empty graph, as expected... (when finer grid of smaller values was used, need enough values close to 0)


# Optimal model found by AIC

ind.opt.aic = which.min(aic.vals)
tau_sq_opt.aic = tau_sq_vals[ind.opt.aic]
mod.opt.aic = fastGHS::fastGHS(x.sf.scaled,tau_sq = tau_sq_opt.aic, epsilon=1e-4,fix_tau=TRUE)

theta.est.aic <- cov2cor(mod.opt.aic$theta)
theta.est.aic[which(abs(theta.est.aic) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.aic!=0)
# 0.007070707
tailoredGlasso::precision(as.matrix(theta.true!=0), theta.est.aic!=0)
# 0.8055556
tailoredGlasso::recall(as.matrix(theta.true!=0), theta.est.aic!=0)
# 0.2929293

tau_sq_opt.aic
# 0.7475374


# A great fit! Worked well. 

# Plot AIC and sparsity

# Skip first val, as the the sparsity got too large (due to small eps => local mode)

df.aic = data.frame(sparsity=sparsities.ebic[-1], AIC = aic.vals[-1], tau_sq=tau_sq_vals[-1])

plot.aic = ggplot2::ggplot(df.aic,  aes(y=AIC,x=tau_sq))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_line(color='deepskyblue2') +geom_point(aes(x=tau_sq_opt.aic,aic.vals[ind.opt.aic]), color='red')

plot.spars = ggplot2::ggplot(df.aic,  aes(y=sparsity,x=tau_sq))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_line(colour='hotpink1')  +geom_point(aes(x=tau_sq_opt.aic,sparsity[ind.opt.aic]), color='red') 

pdf('plots/plots2/AIC_aicscore.pdf')
plot.aic
dev.off()

pdf('plots/plots2/AIC_sparsity.pdf')
plot.spars
dev.off()

# We have marked out the value of tau where the AIC is minimized. 

# EXAMPLE 2 ------------------------------------------------------------------
# GENERATE GRAPH: n=200, p=150

# Test selection by AIC for another scenario

n.2=200
p.2=150
set.seed(12345)
data.sf.2 = huge::huge.generator(n=n.2, d=p.2,graph = 'scale-free') 
g.true.sf.2 = data.sf.2$theta # True adjacency matrix
theta.true.2 = data.sf.2$omega # The precision matrix
theta.true.2[which(theta.true.2<10e-5,arr.ind=T)]=0  
g.sf.2=graph.adjacency(data.sf.2$theta,mode="undirected",diag=F) # true igraph object
x.sf.2 = data.sf.2$data # Observed attributes. nxp matrix.
x.sf.scaled.2= scale(x.sf.2) # Scale columns/variables.
s.sf.scaled.2 = cov(x.sf.scaled.2) # Empirical covariance matrix
data.sf.2$sparsity
# 0.01


# Perform GHS for a range of tau_sq values

tau_sq_vals.2 = seq(1e-5,5,length.out = 40)
sparsities.ebic.2 = rep(0,length(tau_sq_vals.2))
aic.vals.2 = rep(0,length(tau_sq_vals.2))
logliks.ebic.2 = rep(0, length(tau_sq_vals.2))
for(i in 1:length(tau_sq_vals.2)){
  res <- fastGHS::fastGHS(x.sf.scaled.2,tau_sq = tau_sq_vals.2[i],epsilon = 1e-3, fix_tau=TRUE)
  theta.est = res$theta
  theta.est[which(abs(theta.est)<1e-5, arr.ind=T)]=0
  prec.mat = cov2cor(theta.est)
  aic.vals.2[i] = JoStARS::gaussianAIC(s.sf.scaled.2,prec.mat, n.2)
  sparsities.ebic.2[i] = tailoredGlasso::sparsity(theta.est!=0)
  logliks.ebic.2[i] = tailoredGlasso::gaussianloglik(s.sf.scaled.2,prec.mat, n.2)
}


# Optimal model found by AIC

ind.opt.aic.2 = which.min(aic.vals.2)
tau_sq_opt.aic.2 = tau_sq_vals.2[ind.opt.aic.2]
mod.opt.aic.2 = fastGHS::fastGHS(x.sf.scaled.2,tau_sq = tau_sq_opt.aic.2, epsilon=1e-3,fix_tau=TRUE)

theta.est.aic.2 <- cov2cor(mod.opt.aic.2$theta)
theta.est.aic.2[which(abs(theta.est.aic.2) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.aic.2!=0)
# 0.006442953
tailoredGlasso::precision(as.matrix(theta.true.2!=0), theta.est.aic.2!=0)
# 0.8472222
tailoredGlasso::recall(as.matrix(theta.true.2!=0), theta.est.aic.2!=0)
# 0.409396

tau_sq_opt.aic.2
# 4.74359

# Again, good model fit. 

# EXAMPLE 3 ------------------------------------------------------------------
# GENERATE GRAPH: n=200, p=100

# Test selection by AIC for another scenario

n.3=200
p.3=100
set.seed(12345)
data.sf.3 = huge::huge.generator(n=n.3, d=p.3,graph = 'scale-free') 
g.true.sf.3 = data.sf.3$theta # True adjacency matrix
theta.true.3 = data.sf.3$omega # The precision matrix
theta.true.3[which(theta.true.3<10e-5,arr.ind=T)]=0  
g.sf.3=graph.adjacency(data.sf.3$theta,mode="undirected",diag=F) # true igraph object
x.sf.3 = data.sf.3$data # Observed attributes. nxp matrix.
x.sf.scaled.3= scale(x.sf.3) # Scale columns/variables.
s.sf.scaled.3 = cov(x.sf.scaled.3) # Empirical covariance matrix
data.sf.3$sparsity
# 0.01


# Perform GHS for a range of tau_sq values

tau_sq_vals.3 = seq(1e-5,1,length.out = 100)
sparsities.ebic.3 = rep(0,length(tau_sq_vals.3))
aic.vals.3 = rep(0,length(tau_sq_vals.3))
logliks.ebic.3 = rep(0, length(tau_sq_vals.3))
ebic.vals.3 = rep(0,length(tau_sq_vals.3))
for(i in 1:length(tau_sq_vals.3)){
  res <- fastGHS::fastGHS(x.sf.scaled.3,tau_sq = tau_sq_vals.3[i],epsilon = 1e-3, fix_tau=TRUE)
  theta.est = res$theta
  theta.est[which(abs(theta.est)<1e-5, arr.ind=T)]=0
  prec.mat = cov2cor(theta.est)
  aic.vals.3[i] = JoStARS::gaussianAIC(s.sf.scaled.3,prec.mat, n.3)
  sparsities.ebic.3[i] = tailoredGlasso::sparsity(theta.est!=0)
  logliks.ebic.3[i] = tailoredGlasso::gaussianloglik(s.sf.scaled.3,prec.mat, n.3)
  ebic.vals.3[i] = JoStARS::eBIC(s.sf.scaled.3,prec.mat, n.3, gamma=0)
}


# Optimal model found by AIC

ind.opt.aic.3 = which.min(aic.vals.3)
tau_sq_opt.aic.3 = tau_sq_vals.3[ind.opt.aic.3]
mod.opt.aic.3 = fastGHS::fastGHS(x.sf.scaled.3,tau_sq = tau_sq_opt.aic.3, epsilon=1e-3,fix_tau=TRUE)

theta.est.aic.3 <- cov2cor(mod.opt.aic.3$theta)
theta.est.aic.3[which(abs(theta.est.aic.3) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.aic.3!=0)
# 0.005656566
tailoredGlasso::precision(as.matrix(theta.true.3!=0), theta.est.aic.3!=0)
# 0.9642857
tailoredGlasso::recall(as.matrix(theta.true.3!=0), theta.est.aic.3!=0)
# 0.2727273

tau_sq_opt.aic.3
# 0.07071636

# Again, good model fit. 

# EXAMPLE 4 ------------------------------------------------------------------
# GENERATE GRAPH: n=300, p=200

# Test selection by AIC for another scenario

n.4=300
p.4=200
set.seed(12345)
data.sf.4 = huge::huge.generator(n=n.4, d=p.4,graph = 'scale-free') 
g.true.sf.4 = data.sf.4$theta # True adjacency matrix
theta.true.4 = data.sf.4$omega # The precision matrix
theta.true.4[which(theta.true.4<10e-5,arr.ind=T)]=0  
g.sf.4=graph.adjacency(data.sf.4$theta,mode="undirected",diag=F) # true igraph object
x.sf.4 = data.sf.4$data # Observed attributes. nxp matrix.
x.sf.scaled.4= scale(x.sf.4) # Scale columns/variables.
s.sf.scaled.4 = cov(x.sf.scaled.4) # Empirical covariance matrix
data.sf.4$sparsity
# 0.01


# Perform GHS for a range of tau_sq values

tau_sq_vals.4 = seq(1e-5,20,length.out = 80)
sparsities.ebic.4 = rep(0,length(tau_sq_vals.4))
aic.vals.4 = rep(0,length(tau_sq_vals.4))
logliks.ebic.4 = rep(0, length(tau_sq_vals.4))
ebic.vals.4 = rep(0,length(tau_sq_vals.4))
for(i in 1:length(tau_sq_vals.4)){
  res <- fastGHS::fastGHS(x.sf.scaled.4,tau_sq = tau_sq_vals.4[i],epsilon = 1e-3, fix_tau=TRUE)
  theta.est = res$theta
  theta.est[which(abs(theta.est)<1e-5, arr.ind=T)]=0
  prec.mat = cov2cor(theta.est)
  aic.vals.4[i] = JoStARS::gaussianAIC(s.sf.scaled.4,prec.mat, n.4)
  sparsities.ebic.4[i] = tailoredGlasso::sparsity(theta.est!=0)
  logliks.ebic.4[i] = tailoredGlasso::gaussianloglik(s.sf.scaled.4,prec.mat, n.4)
  ebic.vals.4[i] = JoStARS::eBIC(s.sf.scaled.4,prec.mat, n.4, gamma=0)
}


# Optimal model found by AIC

ind.opt.aic.4 = which.min(aic.vals.4)
tau_sq_opt.aic.4 = tau_sq_vals.4[ind.opt.aic.4]
mod.opt.aic.4 = fastGHS::fastGHS(x.sf.scaled.4,tau_sq = tau_sq_opt.aic.4, epsilon=1e-3,fix_tau=TRUE)

theta.est.aic.4 <- cov2cor(mod.opt.aic.4$theta)
theta.est.aic.4[which(abs(theta.est.aic.4) < 1e-5, arr.ind = T)] = 0
tailoredGlasso::sparsity(theta.est.aic.4!=0)
# 0.005678392
tailoredGlasso::precision(as.matrix(theta.true.4!=0), theta.est.aic.4!=0)
# 0.9292035
tailoredGlasso::recall(as.matrix(theta.true.4!=0), theta.est.aic.4!=0)
# 0.5276382

tau_sq_opt.aic.4
# 17.97468

# Again, good model fit. 
