library(jointGHS) # Must be installed from Camiling/jointGHS on github
library(fastGHS) # Must be installed from Camiling/fastGHS
library(tailoredGlasso) # Must be installed from Camiling/tailoredGlasso
library(huge)
library(glasso)
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
source('simulation_functions/perform_jointGHS_simulation.R')
source('simulation_functions/help_functions.R')


# Plot SSL vs jointGHS theta_ij ---------------------------


# n=100, p=100, larger partial correlations (0.203)
n.1=100
p.1=100
set.seed(123)
data.sf.1 = huge::huge.generator(n=n.1, d=p.1,graph = 'scale-free',v=1,u=0.01) 
g.true.sf.1 = data.sf.1$theta # True adjacency matrix
theta.true.1 = data.sf.1$omega # The precision matrix
theta.true.1[which(theta.true.1<10e-5,arr.ind=T)]=0  
x.sf.1 = data.sf.1$data # Observed attributes. nxp matrix.
x.sf.scaled.1= scale(x.sf.1) # Scale columns/variables.
s.sf.scaled.1 = cov(x.sf.scaled.1) # Empirical covariance matrix
data.sf.1$sparsity # True sparsity: 0.04
cov2cor(theta.true.1)[1:5,1:5]          
#[,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000 0.2290578 0.0000000 0.0000000 0.0000000
#[2,] 0.2290578 1.0000000 0.2290578 0.0000000 0.0000000
#[3,] 0.0000000 0.2290578 1.0000000 0.2290578 0.2290578
#[4,] 0.0000000 0.0000000 0.2290578 1.0000000 0.0000000
#[5,] 0.0000000 0.0000000 0.2290578 0.0000000 1.0000000

# Generate second data set 
n.1.2 =  150
set.seed(1233)
graph.1.2 = mutate.graph(data.sf.1, 0.6)
theta.true.1.2 = graph.1.2$prec.mat
x2.sf.1 = mvtnorm::rmvnorm(n.1.2, sigma = cov2cor(graph.1.2$cov.mat))
x2.sf.1.scaled = scale(x2.sf.1)

y=list(x.sf.scaled.1, x2.sf.1.scaled)

penalty <- "fused"
lam1 <- 1
lam2 <- 1
v1 <- 1
lam.eff <- lam1 + c(1:10) * 5
v0s <- lam1/lam.eff
set.seed(123)
fit.ssjgl = SSJGL(Y=y,penalty=penalty,lambda0=1, lambda1=lam1,lambda2=lam2, v1 = v1, v0s = v0s, tol.em=1e-4, a=1, b=p.1, doubly=T, normalize=T)
theta.1.ssjgl = fit.ssjgl$thetalist[[10]][[1]]
theta.2.ssjgl = fit.ssjgl$thetalist[[10]][[2]]
tailoredGlasso::sparsity(theta.1.ssjgl!=0)
# 0.009494949
tailoredGlasso::precision(theta.true.1!=0,theta.1.ssjgl!=0)
# 0.7446809
tailoredGlasso::recall(theta.true.1!=0,theta.1.ssjgl!=0)
# 0.3535354


fit.jointghs = jointGHS::jointGHS(y, AIC_eps=0.1, epsilon = 1e-3)
theta1.est.joint = fit.jointghs$theta[[1]]
theta2.est.joint = fit.jointghs$theta[[2]]
theta1.est.joint[abs(theta1.est.joint)<1e-5]=0
tailoredGlasso::sparsity(theta1.est.joint!=0)
# 0.008686869
tailoredGlasso::precision(theta.true.1!=0,theta1.est.joint!=0)
# 0.9302326
tailoredGlasso::recall(theta.true.1!=0,theta1.est.joint!=0)
# 0.4040404


# Plotting 
theta1.est.joint.upper = cov2cor(theta1.est.joint)[upper.tri(theta1.est.joint)]
theta2.est.joint.upper = cov2cor(theta2.est.joint)[upper.tri(theta2.est.joint)]
theta.1.ssjgl.upper = cov2cor(theta.1.ssjgl)[upper.tri(theta.1.ssjgl)]
theta.2.ssjgl.upper = cov2cor(theta.2.ssjgl)[upper.tri(theta.2.ssjgl)]
truth.all = matrix('none', p.1, p.1)
truth.all[theta.true.1!=0 & theta.true.1.2!=0] = 'both'
truth.all[theta.true.1!=0 & theta.true.1.2==0] = 'network 1'
truth.all[theta.true.1==0 & theta.true.1.2!=0] = 'network 2'
truth.all = truth.all[upper.tri(truth.all)]

df.plot.ssjgl = data.frame(network1=theta.1.ssjgl.upper, network2 = theta.2.ssjgl.upper,truth=factor(truth.all))
df.plot.jointGHS = data.frame(network1=theta1.est.joint.upper, network2=theta2.est.joint.upper,truth=factor(truth.all))
g1 = ggplot2::ggplot(df.plot.jointGHS, aes(y=network2,x=network1))+ labs(title="jointGHS")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=truth, shape=truth)) + geom_hline(yintercept=0, linetype='dashed', color='darkgrey') + geom_vline(xintercept=0,linetype='dashed', color='darkgrey') + theme(legend.position = "none")
g2 = ggplot2::ggplot(df.plot.ssjgl, aes(x =network1, y = network2))+ labs(title="SSJGL")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=truth,shape=truth)) + geom_hline(yintercept=0, linetype='dashed', color='darkgrey') + geom_vline(xintercept=0,linetype='dashed', color='darkgrey')

tmp = ggplot2::ggplot_gtable(ggplot2::ggplot_build(g2))
leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
legend = tmp$grobs[[leg]]
g2 = g2 + theme(legend.position = "none")

pdf('plots/SSJGLvsJointGHS.pdf',13,6)
gridExtra::grid.arrange(g1,g2,legend,ncol=3,widths=c(6,6,1))
dev.off()


# Plot theta_ij vs nu_ij ----------------------------------------------------------

theta.joint = fit.jointghs$theta
Lambda_sq.joint = fit.jointghs$Lambda_sq

# It is important to note that diagonal elements in theta have nothing to do with Lambda. Thus we only consider upper off-diagonal elements (bc of symmetry)
theta.offdiag.1 = abs(cov2cor(theta.joint[[1]]))
theta.offdiag.1 = theta.offdiag.1[upper.tri(theta.offdiag.1)]
theta.offdiag.2 = abs(cov2cor(theta.joint[[2]]))
theta.offdiag.2 = theta.offdiag.2[upper.tri(theta.offdiag.2)]
Lambda_sq.offdiag.1 = Lambda_sq.joint[[1]]
Lambda_sq.offdiag.1 = Lambda_sq.offdiag.1[upper.tri(Lambda_sq.offdiag.1)]
Lambda_sq.offdiag.2 = Lambda_sq.joint[[2]]
Lambda_sq.offdiag.2 = Lambda_sq.offdiag.2[upper.tri(Lambda_sq.offdiag.2)]
# Get E[1/Nu] matrix. Same for all K networks
E_NuInv.iter = fit.jointghs$E_NuInv
E_NuInv.offdiag = E_NuInv.iter[upper.tri(E_NuInv.iter)]

df.plot.2 = data.frame(theta=c(theta.offdiag.1, theta.offdiag.2), Lambda_sq = c(Lambda_sq.offdiag.1, Lambda_sq.offdiag.2), NuInv = rep(E_NuInv.offdiag,2), 
                       truth = factor(rep(truth.all, 2)),estimate= factor(c(rep(1,length(theta.offdiag.1)), rep(2,length(theta.offdiag.2)))))
  
p.lambda = ggplot2::ggplot(df.plot.2, aes(y=Lambda_sq,x=NuInv))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=estimate, shape=truth)) + scale_color_manual(values=c("darkorange", "darkturquoise"))+ 
  geom_hline(yintercept=0, linetype='dashed', color='darkgrey') + geom_vline(xintercept=0,linetype='dashed', color='darkgrey')                     
p.theta = ggplot2::ggplot(df.plot.2, aes(y=theta,x=NuInv))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=estimate, shape=truth))+ scale_color_manual(values=c("darkorange", "darkturquoise")) + 
  geom_hline(yintercept=0, linetype='dashed', color='darkgrey') + geom_vline(xintercept=0,linetype='dashed', color='darkgrey')                

pdf('plots/Theta_vs_NuInvJointGHS.pdf',8,7)
p.theta
dev.off()


