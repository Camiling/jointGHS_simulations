library(fastGHS)
library(jointGHS)
library(huge)
library(glasso)
library(igraph)
library(tailoredGlasso)
library(superheat)
source('simulation_functions/help_functions.R')

# Plot MCC between networks

# First network: n=100, p=100
n.test=100
p.test=50
set.seed(12345)
data.sf.test= huge::huge.generator(n=n.test, d=p.test,graph = 'scale-free',v=0.5,u=0.05) 
g.true.sf.test = data.sf.test$theta # True adjacency matrix
theta.true.test = data.sf.test$omega # The precision matrix
theta.true.test[which(theta.true.test<10e-5,arr.ind=T)]=0  
g.sf.test=graph.adjacency(data.sf.test$theta,mode="undirected",diag=F) # true igraph object
x.sf.test = data.sf.test$data # Observed attributes. nxp matrix.
x.sf.scaled.test= scale(x.sf.test) # Scale columns/variables.
s.sf.scaled.test = cov(x.sf.scaled.test) # Empirical covariance matrix
data.sf.test$sparsity # True sparsity: 0.04

# Generate second data set with same sparsity, 10% edge disagreement
n.test.2=200
p.test.2=50
set.seed(123456)
graph.test.2 = mutate.graph(data.sf.test, 0.1)
theta.true.test.2 = graph.test.2$prec.mat
tailoredGlasso::sparsity(theta.true.test.2!=0)
# 0.04
x.sf.scaled.test.2 = scale(mvtnorm::rmvnorm(n.test.2, sigma = solve(graph.test.2$prec.mat)))

# Generate third data set with same sparsity, 10% edge disagreement
n.test.3=150
p.test.3=50
set.seed(12345677)
graph.test.3 = mutate.graph(data.sf.test, 0.1)
theta.true.test.3 = graph.test.3$prec.mat
tailoredGlasso::sparsity(theta.true.test.3!=0)
# 0.04
x.sf.scaled.test.3 = scale(mvtnorm::rmvnorm(n.test.3, sigma = solve(graph.test.3$prec.mat)))

# Generate fourth data set with same sparsity, unrelated
n.test.4=250
p.test.4=50
set.seed(12345688)
data.sf.test.4 = huge::huge.generator(n=n.test.4, d=p.test.4,graph = 'scale-free',v=0.5,u=0.05) 
g.true.sf.test.4 = data.sf.test.4$theta # True adjacency matrix
theta.true.test.4 = data.sf.test.4$omega # The precision matrix
theta.true.test.4[which(theta.true.test.4<10e-5,arr.ind=T)]=0  
x.sf.test.4 = data.sf.test.4$data # Observed attributes. nxp matrix.
x.sf.scaled.test.4= scale(x.sf.test.4) # Scale columns/variables.
data.sf.test$sparsity # True sparsity: 0.04

thetas.true.list = list(theta.true.test, theta.true.test.2, theta.true.test.3, theta.true.test.4)

# Use jointGHS on the four related data sets
set.seed(1234)
res.joint = jointGHS::jointGHS(list(x.sf.scaled.test, x.sf.scaled.test.2, x.sf.scaled.test.3, x.sf.scaled.test.4), AIC_selection = T, epsilon = 1e-3, AIC_eps = 0.1)

theta1.est <- cov2cor(res.joint$theta[[1]])
theta2.est <- cov2cor(res.joint$theta[[2]])
theta3.est <- cov2cor(res.joint$theta[[3]])
theta4.est <- cov2cor(res.joint$theta[[4]])
theta1.est[which(abs(theta1.est) < 1e-5, arr.ind = T)] = 0
theta2.est[which(abs(theta2.est) < 1e-5, arr.ind = T)] = 0
theta3.est[which(abs(theta3.est) < 1e-5, arr.ind = T)] = 0
theta4.est[which(abs(theta4.est) < 1e-5, arr.ind = T)] = 0

thetas.est.list = list(theta1.est, theta2.est, theta3.est, theta4.est)

# Single version networks

thetas.single.est.list = res.joint$theta_single
thetas.single.est.list = lapply(thetas.single.est.list, FUN = function(s) cov2cor(s))
for(k in 1:length(thetas.single.est.list)){
  thetas.single.est.list[[k]][which(abs(thetas.single.est.list[[k]])<1e-5,arr.ind=T)]=0
}

# Check sparsities
unlist(lapply(thetas.est.list,FUN = function(s) sparsity(s!=0)))

# Plot off-diagonal lambdas
plot(res.joint$Lambda_sq_single[[1]][upper.tri(res.joint$Lambda_sq_single[[1]])], res.joint$Lambda_sq[[1]][upper.tri(res.joint$Lambda_sq[[1]])])

plot(res.joint$Lambda_sq_single[[4]][upper.tri(res.joint$Lambda_sq_single[[4]])], res.joint$Lambda_sq[[4]][upper.tri(res.joint$Lambda_sq[[4]])])

# Plot MCC of joint estimates
MCC_vals = sapply(thetas.est.list, function(x) sapply(thetas.est.list, function(y) MCC(x!=0,y!=0))) # A matrix of pairwisely found MCC
rownames(MCC_vals)=colnames(MCC_vals)= 1:4

png("Figures/MCC_joint_theta.png", width=1200, height =1200)
superheat(MCC_vals,title = 'MCC of jointGHS precision matrix estimates')
dev.off()

# Plot MCC of single estimates 
MCC_vals_single = sapply(thetas.single.est.list, function(x) sapply(thetas.single.est.list, function(y) MCC(x!=0,y!=0))) # A matrix of pairwisely found MCC
rownames(MCC_vals_single)=colnames(MCC_vals_single)= 1:4

png("Figures/MCC_single_theta.png", width=1200, height =1200)
superheat(MCC_vals_single,title = 'MCC of fastGHS precision matrix estimates')
dev.off()

# Plot MCC of single estimates vs joint estimates
MCC_vals_all = sapply(c(thetas.single.est.list, thetas.est.list), function(x) sapply(c(thetas.single.est.list,thetas.est.list), function(y) MCC(x!=0,y!=0))) # A matrix of pairwisely found MCC
rownames(MCC_vals_all)=colnames(MCC_vals_all)= c(paste0(1:4, ' (single)'),paste0(1:4,' (joint)'))

png("Figures/MCC_all_theta.png", width=1200, height =1200)
superheat(MCC_vals_all,title = 'MCC of all precision matrix estimates')
dev.off()


# Plot matrix distance of single-network Lambda_sq
dist_vals = sapply(res.joint$Lambda_sq_single, function(x) sapply(res.joint$Lambda_sq_single, function(y) mean(sqrt((x[upper.tri(x)]-y[upper.tri(y)])^2)))) # A matrix of pairwisely found MCC
rownames(dist_vals )=colnames(dist_vals )= 1:4

png("Figures/Dist_Lambda_sq.png", width=1200, height =1200)
superheat(dist_vals,title = 'Average distance between fastGHS Lambda_sq esimates')
dev.off()


# Plot matrix distance of all Lambda_sq
dist_vals = sapply(c(res.joint$Lambda_sq_single, res.joint$Lambda_sq), function(x) sapply(c(res.joint$Lambda_sq_single,res.joint$Lambda_sq), function(y) mean(sqrt((x[upper.tri(x)]-y[upper.tri(y)])^2)))) # A matrix of pairwisely found MCC
rownames(dist_vals )=colnames(dist_vals )= c(paste0(1:4, ' (single)'),paste0(1:4,' (joint)'))

png("Figures/Dist_Lambda_sq_all.png", width=1200, height =1200)
superheat(dist_vals,title = 'Average distance between all Lambda_sq esimates')
dev.off()
