source('simulation_functions/perform_jointGHS_simulation.R')
source('simulation_functions/help_functions.R')

# Test Williams et al. approach to selecting tau on *denser* network with *larger partial cor* (based on Piironen & Vehtari) for single network -------------------------

# Also check why one network tends to get larger theta values than the other - turns out to be small n that leads to this. 

# Set parameters
p = 50
n=100

set.seed(123)
gg = huge.generator(n=n,d=p, graph='cluster', v=10)
gg$omega[1:10,1:10]
# non-zero off-diag elements around 0.5 in size.
theta.true = (as.matrix(gg$theta)+0)!=0
gg$sparsity
# 0.115102

# Prior expectation of sparsity is the oracle one
alpha.prior = gg$sparsity
tau.prior = alpha.prior/(1-alpha.prior)/sqrt(n)

res.dense = fastGHS::fastGHS(gg$data, tau_sq = tau.prior^2, fix_tau = T)
theta.dense = res.dense$theta
theta.dense[which(abs(theta.dense)<1e-5)] = 0 
sparsity(theta.dense) 
# 0

# Test tau = 1
res.dense.2 = fastGHS::fastGHS(gg$data, tau_sq = 1, fix_tau = T)
theta.dense.2 = res.dense.2$theta
theta.dense.2[which(abs(theta.dense.2)<1e-5)] = 0 
sparsity(theta.dense.2) 
# 0.04326531
precision(theta.true,theta.dense.2!=0)
# 0.7924528

# Compare to glasso
gl = huge(cov(gg$data), lambda=0.787, method='glasso')
gl$sparsity
# 0.04326531
precision(theta.true,gl$icov[[1]]!=0)
# 0.5283019

# Better with GHS with tau fixed. 

# Check why one network gets larger prec matrix elements that the other ------------------------------


nCores=2
N = 2

# Set parameters
K=2
p = 100
n.vals = c(150,200)
stars.thresh = 0.01 # Used in the ordinary graphical lasso
ebic.gamma = 0.2 # Used in JoStARS
var.thresh.jostars = 0.05 # Used in JoStARS

set.seed(1234)
res.inv.1 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = 0.2, tau_sq = 10, tau_sq_ghs = 0.01,  ebic.gamma=ebic.gamma,stars.thresh=stars.thresh, 
                                        var.thresh.jostars = var.thresh.jostars,method='symmetric', nCores = nCores, save.res.jointGHS = T)


theta.true.1.1 = cov2cor(res.inv.1$true.prec.matrices[[1]])
theta.true.1.1 = theta.true.1.1[upper.tri(theta.true.1.1)]
theta.true.1.2 = res.inv.1$true.prec.matrices[[2]]
theta.true.1.2 = theta.true.1.2[upper.tri(theta.true.1.2)]

dev.new()
par(mfrow=c(1,2))
hist(theta.true.1.1, breaks=100)
hist(theta.true.1.2, breaks=100)

# First prec mat is unscaled => must scale before comparing. Then they are identical. 

# Why are the elements much larger in estimate?

p=100
huge.init = huge.generator(n,p,graph='scale-free',verbose = F,v=0.5,u=0.05)
#theta.mut = mutate.graph(huge.init,0.2,scale=T)

x1 = scale(mvtnorm::rmvnorm(n, mean=rep(0,p), huge.init$sigma))
x2 = mvtnorm::rmvnorm(n, mean=rep(0,p), cov2cor(huge.init$sigma))

s1 = cov(x1)
s2 = cov(x2)
s1=s1[upper.tri(s1)]
s2=s2[upper.tri(s2)]

dev.new()
par(mfrow=c(1,2))
hist(s1, breaks=100)
hist(s2, breaks=100)

# Seems identical - cov2cor before or scale after should give same res....

# Test other seed

p=100
set.seed(12)
res.inv.2 = perform_jointGHS_simulation(K,n.vals=c(200,150), p, N, frac.disagreement = 0.2, tau_sq = 10, tau_sq_ghs = 0.01,  ebic.gamma=ebic.gamma,stars.thresh=stars.thresh, 
                                        var.thresh.jostars = var.thresh.jostars,method='symmetric', nCores = nCores, save.res.jointGHS = T)


# look at non-zero elements of theta
theta.2.1 = res.inv.2$theta[[1]][[1]]
theta.2.1 = theta.2.1[upper.tri(theta.2.1)]
theta.2.1 = theta.2.1[which(abs(theta.2.1)>1e-5)]
theta.2.2 = res.inv.2$theta[[1]][[2]]
theta.2.2 = theta.2.2[upper.tri(theta.2.2)]
theta.2.2 = theta.2.2[which(abs(theta.2.2)>1e-5)]

dev.new()
par(mfrow=c(1,2))
hist(theta.2.1, breaks=100)
hist(theta.2.2, breaks=100)


# Create plot of theta vs lambda
theta.iter.2.1 = res.inv.2$theta[[1]]
Lambda_sq.iter.2.1 = res.inv.2$Lambda_sq[[1]]

# It is important to note that diagonal elements in theta have nothing to do with Lambda. Thus we only consider upper off-diagonal elements (bc of symmetry)
theta.offdiag.2.1.1 = theta.iter.2.1[[1]]
theta.offdiag.2.1.1[!upper.tri(theta.offdiag.2.1.1)] = NA
theta.offdiag.2.1.2 = theta.iter.2.1[[2]]
theta.offdiag.2.1.2[!upper.tri(theta.offdiag.2.1.2)] = NA
Lambda_sq.offdiag.2.1.1 = Lambda_sq.iter.2.1[[1]]
Lambda_sq.offdiag.2.1.1[!upper.tri(Lambda_sq.offdiag.2.1.1)] = NA
Lambda_sq.offdiag.2.1.2 = Lambda_sq.iter.2.1[[2]]
Lambda_sq.offdiag.2.1.2[!upper.tri(Lambda_sq.offdiag.2.1.2)] = NA
# Get E[1/Nu] matrix. Same for all K networks
E_NuInv.iter.2.1 = res.inv.2$E_NuInv[[1]]
E_NuInv.offdiag.2.1 = E_NuInv.iter.2.1
E_NuInv.offdiag.2.1[!upper.tri(E_NuInv.offdiag.2.1)] = NA

# Add information on whether an edge is present in network 1, network 2 or both.
mark.mat = matrix('none',p, p)
mark.mat[res.inv.2$true.prec.matrices[[1]]!=0] = '1'
mark.mat[res.inv.2$true.prec.matrices[[2]] !=0] = '2'
mark.mat[res.inv.2$true.prec.matrices[[1]]!=0 & res.inv.2$true.prec.matrices[[2]] !=0] = 'both'

df.all.iter.2 = data.frame(theta=c(c(theta.offdiag.2.1.1), c(theta.offdiag.2.1.2)), Lambda = c(c(sqrt(Lambda_sq.offdiag.2.1.1)), c(sqrt(Lambda_sq.offdiag.2.1.2))), 
                           NuInv=c(c(E_NuInv.offdiag.2.1), c(E_NuInv.offdiag.2.1)),
                           truth=factor(c(c(mark.mat), c(mark.mat))),
                           graph = factor(c(rep(1,p^2), rep(2, p^2))))


# Plot theta_ijk's as functions of lambda_ijk's

p1.2 <- ggplot2::ggplot(na.omit(df.all.iter.2), aes(y=theta,x=Lambda))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=graph, shape=truth))     

p1.2

# Now the opposite - seems like smaller n => larger theta vals

