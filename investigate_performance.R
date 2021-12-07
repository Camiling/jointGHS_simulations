source('simulation_functions/perform_jointGHS_simulation.R')
source('simulation_functions/help_functions.R')
library(ggplot2)

# Investigate to behaviour of the jointGHS when fixing tau to a sensible value
# Create Lambda vs NuInv plots, theta vs Lambda plots etc. 


nCores=2
N = 2

# Set parameters
K=2
p = 100
n.vals = c(150,200)
stars.thresh = 0.01 # Used in the ordinary graphical lasso
ebic.gamma = 0.2 # Used in JoStARS
var.thresh.jostars = 0.05 # Used in JoStARS

# Case 1: K=2, similar distributions (80% edge agreement) --------------------

# Datasets from similar distributions (80% edge agreement)

set.seed(1234)
res.inv.1 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = 0.2, tau_sq = 10, tau_sq_ghs = 0.01,  ebic.gamma=ebic.gamma,stars.thresh=stars.thresh, 
                                      var.thresh.jostars = var.thresh.jostars,method='symmetric', nCores = nCores, save.res.jointGHS = T)

# For first case, first replicate 
theta.iter.1.1 = res.inv.1$theta[[1]]
Lambda_sq.iter.1.1 = res.inv.1$Lambda_sq[[1]]

# It is important to note that diagonal elements in theta have nothing to do with Lambda. Thus we only consider upper off-diagonal elements (bc of symmetry)
theta.offdiag.1.1.1 = theta.iter.1.1[[1]]
theta.offdiag.1.1.1[!upper.tri(theta.offdiag.1.1.1)] = NA
theta.offdiag.1.1.2 = theta.iter.1.1[[2]]
theta.offdiag.1.1.2[!upper.tri(theta.offdiag.1.1.2)] = NA
Lambda_sq.offdiag.1.1.1 = Lambda_sq.iter.1.1[[1]]
Lambda_sq.offdiag.1.1.1[!upper.tri(Lambda_sq.offdiag.1.1.1)] = NA
Lambda_sq.offdiag.1.1.2 = Lambda_sq.iter.1.1[[2]]
Lambda_sq.offdiag.1.1.2[!upper.tri(Lambda_sq.offdiag.1.1.2)] = NA
# Get E[1/Nu] matrix. Same for all K networks
E_NuInv.iter.1.1 = res.inv.1$E_NuInv[[1]]
E_NuInv.offdiag.1.1 = E_NuInv.iter.1.1
E_NuInv.offdiag.1.1[!upper.tri(E_NuInv.offdiag.1.1)] = NA

# Add information on whether an edge is present in network 1, network 2 or both.
mark.mat = matrix('none',p, p)
mark.mat[res.inv.1$true.prec.matrices[[1]]!=0] = '1'
mark.mat[res.inv.1$true.prec.matrices[[2]] !=0] = '2'
mark.mat[res.inv.1$true.prec.matrices[[1]]!=0 & res.inv.1$true.prec.matrices[[2]] !=0] = 'both'

df.all.iter.1 = data.frame(theta=c(c(theta.offdiag.1.1.1), c(theta.offdiag.1.1.2)), Lambda = c(c(sqrt(Lambda_sq.offdiag.1.1.1)), c(sqrt(Lambda_sq.offdiag.1.1.2))), 
                           NuInv=c(c(E_NuInv.offdiag.1.1), c(E_NuInv.offdiag.1.1)),
                           truth=factor(c(c(mark.mat), c(mark.mat))),
                           graph = factor(c(rep(1,p^2), rep(2, p^2))))


# Plot theta_ijk's as functions of lambda_ijk's


p1 <- ggplot2::ggplot(na.omit(df.all.iter.1), aes(y=theta,x=Lambda))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
                   geom_point(aes(colour=graph))                      

pdf('plots/theta_vs_Lambda_K2_80sim.pdf')
p1
dev.off()

p1.truth <- ggplot2::ggplot(na.omit(df.all.iter.1), aes(y=theta,x=Lambda))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=graph, shape=truth))                      

pdf('plots/theta_vs_Lambda_K2_80sim_truth.pdf')
p1.truth
dev.off()


# Plot lambda_ijk's as functions of nu_ij's

# Using the same scenario

p2 <- ggplot2::ggplot(na.omit(df.all.iter.1), aes(y=Lambda,x=NuInv))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=graph))                      

pdf('plots/Lambda_vs_NuInv_K2_80sim.pdf')
p2
dev.off()

p2.truth <- ggplot2::ggplot(na.omit(df.all.iter.1), aes(y=Lambda,x=NuInv))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=graph, shape=truth))                      

pdf('plots/Lambda_vs_NuInv_K2_80sim_truth.pdf')
p2.truth
dev.off()

# Plot theta_ijk's as functions of nu_ij's

p3 <- ggplot2::ggplot(na.omit(df.all.iter.1), aes(y=theta,x=NuInv))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=graph))                      

pdf('plots/theta_vs_NuInv_K2_80sim.pdf')
p3
dev.off()

p3.truth <- ggplot2::ggplot(na.omit(df.all.iter.1), aes(y=theta,x=NuInv))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=graph, shape=truth))                      

pdf('plots/theta_vs_NuInv_K2_80sim_truth.pdf')
p3.truth
dev.off()



# Case 2: K=2, identical distributions (100% edge agreement) ----------------------

set.seed(1234)
res.inv.2 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = 0, tau_sq = 10, tau_sq_ghs = 0.01,  ebic.gamma=ebic.gamma,stars.thresh=stars.thresh, 
                                        var.thresh.jostars = var.thresh.jostars,method='symmetric', nCores = nCores, save.res.jointGHS = T)

# For second case, first replicate 
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
                           truth = factor(c(c(mark.mat), c(mark.mat))), 
                           graph = factor(c(rep(1,p^2), rep(2, p^2))))


# Plot theta_ijk's as functions of lambda_ijk's


p1.2 <- ggplot2::ggplot(na.omit(df.all.iter.2), aes(y=theta,x=Lambda))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=graph))                      

pdf('plots/theta_vs_Lambda_K2_100sim.pdf')
p1.2
dev.off()

p1.2.truth <- ggplot2::ggplot(na.omit(df.all.iter.2),  aes(y=theta,x=Lambda))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=graph, shape=truth))                      

pdf('plots/theta_vs_Lambda_K2_100sim_truth.pdf')
p1.2.truth 
dev.off()


# Plot lambda_ijk's as functions of nu_ij's

# Using the same scenario

p2.2 <- ggplot2::ggplot(na.omit(df.all.iter.2), aes(y=Lambda,x=NuInv))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=graph))                      

pdf('plots/Lambda_vs_NuInv_K2_100sim.pdf')
p2.2
dev.off()

p2.2.truth <- ggplot2::ggplot(na.omit(df.all.iter.2),  aes(y=Lambda,x=NuInv))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=graph, shape=truth))                      

pdf('plots/Lambda_vs_NuInv_K2_100sim_truth.pdf')
p2.2.truth 
dev.off()

# Plot theta_ijk's as functions of nu_ij's

p3.2 <- ggplot2::ggplot(na.omit(df.all.iter.2), aes(y=theta,x=NuInv))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=graph))                      

pdf('plots/theta_vs_NuInv_K2_100sim.pdf')
p3.2
dev.off()


p3.2.truth <- ggplot2::ggplot(na.omit(df.all.iter.2),  aes(y=theta,x=NuInv))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=graph, shape=truth))                      

pdf('plots/theta_vs_NuInv_K2_100sim_truth.pdf')
p3.2.truth 
dev.off()





# Case 3: K=2, unrelated distributions (0% edge agreement) ----------------------

set.seed(123)
res.inv.3 = perform_jointGHS_simulation(K,n.vals, p, N, frac.disagreement = 1, tau_sq = 10, tau_sq_ghs = 0.01,  ebic.gamma=ebic.gamma,stars.thresh=stars.thresh, 
                                        var.thresh.jostars = var.thresh.jostars,method='symmetric', nCores = nCores, save.res.jointGHS = T)

# For second case, first replicate 
theta.iter.3.1 = res.inv.3$theta[[1]]
Lambda_sq.iter.3.1 = res.inv.3$Lambda_sq[[1]]

# It is important to note that diagonal elements in theta have nothing to do with Lambda. Thus we only consider upper off-diagonal elements (bc of symmetry)
theta.offdiag.3.1.1 = theta.iter.3.1[[1]]
theta.offdiag.3.1.1[!upper.tri(theta.offdiag.3.1.1)] = NA
theta.offdiag.3.1.2 = theta.iter.3.1[[2]]
theta.offdiag.3.1.2[!upper.tri(theta.offdiag.3.1.2)] = NA
Lambda_sq.offdiag.3.1.1 = Lambda_sq.iter.3.1[[1]]
Lambda_sq.offdiag.3.1.1[!upper.tri(Lambda_sq.offdiag.3.1.1)] = NA
Lambda_sq.offdiag.3.1.2 = Lambda_sq.iter.3.1[[2]]
Lambda_sq.offdiag.3.1.2[!upper.tri(Lambda_sq.offdiag.3.1.2)] = NA
# Get E[1/Nu] matrix. Same for all K networks
E_NuInv.iter.3.1 = res.inv.3$E_NuInv[[1]]
E_NuInv.offdiag.3.1 = E_NuInv.iter.3.1
E_NuInv.offdiag.3.1[!upper.tri(E_NuInv.offdiag.3.1)] = NA

# Add information on whether an edge is present in network 1, network 2 or both.
mark.mat = matrix('none',p, p)
mark.mat[res.inv.3$true.prec.matrices[[1]]!=0] = '1'
mark.mat[res.inv.3$true.prec.matrices[[2]] !=0] = '2'
mark.mat[res.inv.3$true.prec.matrices[[1]]!=0 & res.inv.3$true.prec.matrices[[2]] !=0] = 'both'

df.all.iter.3 = data.frame(theta=c(c(theta.offdiag.3.1.1), c(theta.offdiag.3.1.2)), Lambda = c(c(sqrt(Lambda_sq.offdiag.3.1.1)), c(sqrt(Lambda_sq.offdiag.3.1.2))), 
                           NuInv=c(c(E_NuInv.offdiag.3.1), c(E_NuInv.offdiag.3.1)),
                           truth = factor(c(c(mark.mat), c(mark.mat))), 
                           graph = factor(c(rep(1,p^2), rep(2, p^2))))


# Plot theta_ijk's as functions of lambda_ijk's


p1.3 <- ggplot2::ggplot(na.omit(df.all.iter.3), aes(y=theta,x=Lambda))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=graph))                      

pdf('plots/theta_vs_Lambda_K2_0sim.pdf')
p1.3
dev.off()

p1.3.truth <- ggplot2::ggplot(na.omit(df.all.iter.3),  aes(y=theta,x=Lambda))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=graph, shape=truth))                      

pdf('plots/theta_vs_Lambda_K2_0sim_truth.pdf')
p1.3.truth 
dev.off()



# Plot lambda_ijk's as functions of nu_ij's

# Using the same scenario

p2.3 <- ggplot2::ggplot(na.omit(df.all.iter.3), aes(y=Lambda,x=NuInv))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=graph))                      

pdf('plots/Lambda_vs_NuInv_K2_0sim.pdf')
p2.3
dev.off()

#REDO this plot
p2.3.truth <- ggplot2::ggplot(na.omit(df.all.iter.3),  aes(y=Lambda,x=NuInv))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=graph, shape=truth))                      

pdf('plots/Lambda_vs_NuInv_K2_0sim_truth.pdf')
p2.3.truth 
dev.off()

# Plot theta_ijk's as functions of nu_ij's

p3.3 <- ggplot2::ggplot(na.omit(df.all.iter.3), aes(y=theta,x=NuInv))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=graph))                      

pdf('plots/theta_vs_NuInv_K2_0sim.pdf')
p3.3
dev.off()

p3.3.truth <- ggplot2::ggplot(na.omit(df.all.iter.3),  aes(y=theta,x=Lambda))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=graph, shape=truth))                      

pdf('plots/theta_vs_NuInv_K2_0sim_truth.pdf')
p3.3.truth 
dev.off()





# Case 4: K=4, identical distributions (80% edge agreement), one out -----------------

# Datasets from similar distributions (80% edge agreement)

K.4= 4
n.vals.4.out = c(150, 200, 150, 200)
set.seed(1234)
res.inv.4 = perform_jointGHS_simulation(K.4,n.vals.4.out, p, N, frac.disagreement = 0.2, tau_sq = 5, tau_sq_ghs = 0.01,  ebic.gamma=ebic.gamma,stars.thresh=stars.thresh, 
                                        include.jostars = F, include.JGL = F,method='unsymmetric', nCores = nCores, save.res.jointGHS = T)

# For first case, first replicate 
theta.iter.4.1 = res.inv.4$theta[[1]]
Lambda_sq.iter.4.1 = res.inv.4$Lambda_sq[[1]]

# It is important to note that diagonal elements in theta have nothing to do with Lambda. Thus we only consider upper off-diagonal elements (bc of symmetry)
theta.offdiag.4.1.1 = theta.iter.4.1[[1]]
theta.offdiag.4.1.1[!upper.tri(theta.offdiag.4.1.1)] = NA
theta.offdiag.4.1.2 = theta.iter.4.1[[2]]
theta.offdiag.4.1.2[!upper.tri(theta.offdiag.4.1.2)] = NA
theta.offdiag.4.1.3 = theta.iter.4.1[[3]]
theta.offdiag.4.1.3[!upper.tri(theta.offdiag.4.1.3)] = NA
theta.offdiag.4.1.4 = theta.iter.4.1[[4]]
theta.offdiag.4.1.4[!upper.tri(theta.offdiag.4.1.4)] = NA
Lambda_sq.offdiag.4.1.1 = Lambda_sq.iter.4.1[[1]]
Lambda_sq.offdiag.4.1.1[!upper.tri(Lambda_sq.offdiag.4.1.1)] = NA
Lambda_sq.offdiag.4.1.2 = Lambda_sq.iter.4.1[[2]]
Lambda_sq.offdiag.4.1.2[!upper.tri(Lambda_sq.offdiag.4.1.2)] = NA
Lambda_sq.offdiag.4.1.3 = Lambda_sq.iter.4.1[[3]]
Lambda_sq.offdiag.4.1.3[!upper.tri(Lambda_sq.offdiag.4.1.3)] = NA
Lambda_sq.offdiag.4.1.4 = Lambda_sq.iter.4.1[[4]]
Lambda_sq.offdiag.4.1.4[!upper.tri(Lambda_sq.offdiag.4.1.4)] = NA
# Get E[1/Nu] matrix. Same for all K networks
E_NuInv.iter.4.1 = res.inv.4$E_NuInv[[1]]
E_NuInv.offdiag.4.1 = E_NuInv.iter.4.1
E_NuInv.offdiag.4.1[!upper.tri(E_NuInv.offdiag.4.1)] = NA

# Add information on whether an edge is present in all networks or not. 
mark.mat = matrix('none',p, p)
sum.mat = (res.inv.4$true.prec.matrices[[1]]!=0) + (res.inv.4$true.prec.matrices[[2]]!=0) + (res.inv.4$true.prec.matrices[[3]]!=0) + (res.inv.4$true.prec.matrices[[4]]!=0)
mark.mat[sum.mat==1] = 'one'
mark.mat[sum.mat==2] = 'two'
mark.mat[sum.mat==3] = 'three'
mark.mat[sum.mat==4] = 'all'

df.all.iter.4 = data.frame(theta=c(c(theta.offdiag.4.1.1), c(theta.offdiag.4.1.2), c(theta.offdiag.4.1.3), c(theta.offdiag.4.1.4)), 
                           Lambda = c(c(sqrt(Lambda_sq.offdiag.4.1.1)), c(sqrt(Lambda_sq.offdiag.4.1.2)), c(sqrt(Lambda_sq.offdiag.4.1.3)), c(sqrt(Lambda_sq.offdiag.4.1.4))), 
                           NuInv=c(c(E_NuInv.offdiag.4.1), c(E_NuInv.offdiag.4.1),  c(E_NuInv.offdiag.4.1),  c(E_NuInv.offdiag.4.1)),
                           truth = factor(c(c(mark.mat), c(mark.mat), c(mark.mat), c(mark.mat))), 
                           graph = factor(c(rep(1,p^2), rep(2, p^2), rep(3, p^2), rep(4, p^2))))


# Plot theta_ijk's as functions of lambda_ijk's


p1.4 <- ggplot2::ggplot(na.omit(df.all.iter.4), aes(y=theta,x=Lambda))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=graph))                      

pdf('plots/theta_vs_Lambda_K4_oneeout_80sim.pdf')
p1.4
dev.off()

p1.4.truth <- ggplot2::ggplot(na.omit(df.all.iter.4),  aes(y=theta,x=Lambda))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=graph, shape=truth))                      

pdf('plots/theta_vs_Lambda_K4_oneeout_80sim_truth.pdf')
p1.4.truth 
dev.off()


# Plot lambda_ijk's as functions of nu_ij's

# Using the same scenario

p2.4 <- ggplot2::ggplot(na.omit(df.all.iter.4), aes(y=Lambda,x=NuInv))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=graph))                      

pdf('plots/Lambda_vs_NuInv_K4_oneeout_80sim.pdf')
p2.4
dev.off()

p2.4.truth <- ggplot2::ggplot(na.omit(df.all.iter.4),  aes(y=Lambda,x=NuInv))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=graph, shape=truth))                      

pdf('plots/Lambda_vs_NuInv_K4_oneeout_80sim_truth.pdf')
p2.4.truth 
dev.off()

# Plot theta_ijk's as functions of nu_ij's

p3.4 <- ggplot2::ggplot(na.omit(df.all.iter.4), aes(y=theta,x=NuInv))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=graph))                      

pdf('plots/theta_vs_NuInv_K4_oneeout_80sim.pdf')
p3.4
dev.off()

p3.4.truth <- ggplot2::ggplot(na.omit(df.all.iter.4),  aes(y=theta,x=NuInv))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=graph, shape=truth))                      

pdf('plots/theta_vs_NuInv_K4_oneeout_80sim_truth.pdf')
p3.4.truth 
dev.off()

# REMOVE: 
# Closer look at the negative theta values in graphs 1 and 4
#edge.neg.1 = which(theta.iter.4.1[[1]]< -(1e-5), arr.ind=T) 
#edge.neg.1
# Check which ones it is present in
#lapply(res.inv.4$true.prec.matrices, FUN = function(m) (abs(m[edge.neg.1])>1e-5))
#res.inv.4$true.prec.matrices[[1]][edge.neg.1]
#res.inv.4$true.prec.matrices[[2]][edge.neg.1]
#res.inv.4$true.prec.matrices[[3]][edge.neg.1]
#res.inv.4$true.prec.matrices[[4]][edge.neg.1]
#res.inv.4$E_NuInv[[1]][53, 17]
#edge.neg.4 = which(theta.iter.4.1[[4]]< -(1e-5), arr.ind=T)
#edge.neg.4


# Case 5: K=4, identical distributions (100% edge agreement), one out -----------------

# Datasets from similar distributions (100% edge agreement)

K.5= 4
n.vals.5.out = c(150, 150, 150, 150)
set.seed(1234)
res.inv.5 = perform_jointGHS_simulation(K.5,n.vals.5.out, p, N, frac.disagreement = 0, tau_sq = 5, tau_sq_ghs = 0.01,  ebic.gamma=ebic.gamma,stars.thresh=stars.thresh, 
                                        include.jostars = F, include.JGL = F, method='unsymmetric', nCores = nCores, save.res.jointGHS = T)

# For first case, first replicate 
theta.iter.5.1 = res.inv.5$theta[[1]]
Lambda_sq.iter.5.1 = res.inv.5$Lambda_sq[[1]]

# It is important to note that diagonal elements in theta have nothing to do with Lambda. Thus we only consider upper off-diagonal elements (bc of symmetry)
theta.offdiag.5.1.1 = theta.iter.5.1[[1]]
theta.offdiag.5.1.1[!upper.tri(theta.offdiag.5.1.1)] = NA
theta.offdiag.5.1.2 = theta.iter.5.1[[2]]
theta.offdiag.5.1.2[!upper.tri(theta.offdiag.5.1.2)] = NA
theta.offdiag.5.1.3 = theta.iter.5.1[[3]]
theta.offdiag.5.1.3[!upper.tri(theta.offdiag.5.1.3)] = NA
theta.offdiag.5.1.4 = theta.iter.5.1[[4]]
theta.offdiag.5.1.4[!upper.tri(theta.offdiag.5.1.4)] = NA
Lambda_sq.offdiag.5.1.1 = Lambda_sq.iter.5.1[[1]]
Lambda_sq.offdiag.5.1.1[!upper.tri(Lambda_sq.offdiag.5.1.1)] = NA
Lambda_sq.offdiag.5.1.2 = Lambda_sq.iter.5.1[[2]]
Lambda_sq.offdiag.5.1.2[!upper.tri(Lambda_sq.offdiag.5.1.2)] = NA
Lambda_sq.offdiag.5.1.3 = Lambda_sq.iter.5.1[[3]]
Lambda_sq.offdiag.5.1.3[!upper.tri(Lambda_sq.offdiag.5.1.3)] = NA
Lambda_sq.offdiag.5.1.4 = Lambda_sq.iter.5.1[[4]]
Lambda_sq.offdiag.5.1.4[!upper.tri(Lambda_sq.offdiag.5.1.4)] = NA
# Get E[1/Nu] matrix. Same for all K networks
E_NuInv.iter.5.1 = res.inv.5$E_NuInv[[1]]
E_NuInv.offdiag.5.1 = E_NuInv.iter.5.1
E_NuInv.offdiag.5.1[!upper.tri(E_NuInv.offdiag.5.1)] = NA

# Add information on whether an edge is present in all networks or not. 
mark.mat = matrix('none',p, p)
sum.mat = (res.inv.5$true.prec.matrices[[1]]!=0) + (res.inv.5$true.prec.matrices[[2]]!=0) + (res.inv.5$true.prec.matrices[[3]]!=0) + (res.inv.5$true.prec.matrices[[4]]!=0)
mark.mat[sum.mat==1] = 'one'
mark.mat[sum.mat==2] = 'two'
mark.mat[sum.mat==3] = 'three'
mark.mat[sum.mat==4] = 'all'

df.all.iter.5 = data.frame(theta=c(c(theta.offdiag.5.1.1), c(theta.offdiag.5.1.2), c(theta.offdiag.5.1.3), c(theta.offdiag.5.1.4)), 
                           Lambda = c(c(sqrt(Lambda_sq.offdiag.5.1.1)), c(sqrt(Lambda_sq.offdiag.5.1.2)), c(sqrt(Lambda_sq.offdiag.5.1.3)), c(sqrt(Lambda_sq.offdiag.5.1.4))), 
                           NuInv=c(c(E_NuInv.offdiag.5.1), c(E_NuInv.offdiag.5.1),  c(E_NuInv.offdiag.5.1),  c(E_NuInv.offdiag.5.1)),
                           truth = factor(c(c(mark.mat), c(mark.mat), c(mark.mat), c(mark.mat))), 
                           graph = factor(c(rep(1,p^2), rep(2, p^2), rep(3, p^2), rep(4, p^2))))


# Plot theta_ijk's as functions of lambda_ijk's


p1.5 <- ggplot2::ggplot(na.omit(df.all.iter.5), aes(y=theta,x=Lambda))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=graph))                      

pdf('plots/theta_vs_Lambda_K4_oneeout_100sim.pdf')
p1.5
dev.off()

p1.5.truth <- ggplot2::ggplot(na.omit(df.all.iter.5),  aes(y=theta,x=Lambda))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=graph, shape=truth))                      

pdf('plots/theta_vs_Lambda_K4_oneeout_100sim_truth.pdf')
p1.5.truth 
dev.off()


# Plot lambda_ijk's as functions of nu_ij's

# Using the same scenario

p2.5 <- ggplot2::ggplot(na.omit(df.all.iter.5), aes(y=Lambda,x=NuInv))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=graph))                      

pdf('plots/Lambda_vs_NuInv_K4_oneeout_100sim.pdf')
p2.5
dev.off()

p2.4.truth <- ggplot2::ggplot(na.omit(df.all.iter.5),  aes(y=Lambda,x=NuInv))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=graph, shape=truth))                      

pdf('plots/Lambda_vs_NuInv_K4_oneeout_100sim_truth.pdf')
p2.4.truth 
dev.off()

# Plot theta_ijk's as functions of nu_ij's

p3.5 <- ggplot2::ggplot(na.omit(df.all.iter.5), aes(y=theta,x=NuInv))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=graph))                      

pdf('plots/theta_vs_NuInv_K4_oneeout_100sim.pdf')
p3.5
dev.off()

p3.5.truth <- ggplot2::ggplot(na.omit(df.all.iter.5),  aes(y=theta,x=NuInv))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=graph, shape=truth))                      

pdf('plots/theta_vs_NuInv_K4_oneeout_100sim_truth.pdf')
p3.5.truth 
dev.off()



# Case 6: K=6, identical distributions (100% edge agreement), one out -----------------

# Datasets from similar distributions (100% edge agreement)

K.6= 6
n.vals.6.out = c(150, 150, 150, 150, 150, 150)
set.seed(12345)
res.inv.6 = perform_jointGHS_simulation(K.6,n.vals.6.out, p, N, frac.disagreement = 0, tau_sq = 10, tau_sq_ghs = 0.01,  ebic.gamma=ebic.gamma,stars.thresh=stars.thresh, 
                                        include.jostars = F, include.JGL = F, method='unsymmetric', nCores = nCores, save.res.jointGHS = T)

# For first case, first replicate (N=1)
theta.iter.6.1 = res.inv.6$theta[[1]]
Lambda_sq.iter.6.1 = res.inv.6$Lambda_sq[[1]]

# It is important to note that diagonal elements in theta have nothing to do with Lambda. Thus we only consider upper off-diagonal elements (bc of symmetry)
theta.offdiag.6.1.1 = theta.iter.6.1[[1]]
theta.offdiag.6.1.1[!upper.tri(theta.offdiag.6.1.1)] = NA
theta.offdiag.6.1.2 = theta.iter.6.1[[2]]
theta.offdiag.6.1.2[!upper.tri(theta.offdiag.6.1.2)] = NA
theta.offdiag.6.1.3 = theta.iter.6.1[[3]]
theta.offdiag.6.1.3[!upper.tri(theta.offdiag.6.1.3)] = NA
theta.offdiag.6.1.4 = theta.iter.6.1[[4]]
theta.offdiag.6.1.4[!upper.tri(theta.offdiag.6.1.4)] = NA
theta.offdiag.6.1.5 = theta.iter.6.1[[5]]
theta.offdiag.6.1.5[!upper.tri(theta.offdiag.6.1.5)] = NA
theta.offdiag.6.1.6 = theta.iter.6.1[[6]]
theta.offdiag.6.1.6[!upper.tri(theta.offdiag.6.1.6)] = NA
Lambda_sq.offdiag.6.1.1 = Lambda_sq.iter.6.1[[1]]
Lambda_sq.offdiag.6.1.1[!upper.tri(Lambda_sq.offdiag.6.1.1)] = NA
Lambda_sq.offdiag.6.1.2 = Lambda_sq.iter.6.1[[2]]
Lambda_sq.offdiag.6.1.2[!upper.tri(Lambda_sq.offdiag.6.1.2)] = NA
Lambda_sq.offdiag.6.1.3 = Lambda_sq.iter.6.1[[3]]
Lambda_sq.offdiag.6.1.3[!upper.tri(Lambda_sq.offdiag.6.1.3)] = NA
Lambda_sq.offdiag.6.1.4 = Lambda_sq.iter.6.1[[4]]
Lambda_sq.offdiag.6.1.4[!upper.tri(Lambda_sq.offdiag.6.1.4)] = NA
Lambda_sq.offdiag.6.1.5 = Lambda_sq.iter.6.1[[5]]
Lambda_sq.offdiag.6.1.5[!upper.tri(Lambda_sq.offdiag.6.1.5)] = NA
Lambda_sq.offdiag.6.1.6 = Lambda_sq.iter.6.1[[6]]
Lambda_sq.offdiag.6.1.6[!upper.tri(Lambda_sq.offdiag.6.1.6)] = NA
# Get E[1/Nu] matrix. Same for all K networks
E_NuInv.iter.6.1 = res.inv.6$E_NuInv[[1]]
E_NuInv.offdiag.6.1 = E_NuInv.iter.6.1
E_NuInv.offdiag.6.1[!upper.tri(E_NuInv.offdiag.6.1)] = NA

# Add information on whether an edge is present in all networks or not. 
mark.mat = matrix('none',p, p)
sum.mat = (res.inv.6$true.prec.matrices[[1]]!=0) + (res.inv.6$true.prec.matrices[[2]]!=0) + (res.inv.6$true.prec.matrices[[3]]!=0) + 
    (res.inv.6$true.prec.matrices[[4]]!=0) + (res.inv.6$true.prec.matrices[[5]]!=0) + (res.inv.6$true.prec.matrices[[6]]!=0)
mark.mat[sum.mat==1] = 'one'
mark.mat[sum.mat==2] = 'two'
mark.mat[sum.mat==3] = 'three'
mark.mat[sum.mat==4] = 'four'
mark.mat[sum.mat==5] = 'five'
mark.mat[sum.mat==6] = 'all'

df.all.iter.6 = data.frame(theta=c(c(theta.offdiag.6.1.1), c(theta.offdiag.6.1.2), c(theta.offdiag.6.1.3), c(theta.offdiag.6.1.4), c(theta.offdiag.6.1.5), c(theta.offdiag.6.1.6)), 
                           Lambda = c(c(sqrt(Lambda_sq.offdiag.6.1.1)), c(sqrt(Lambda_sq.offdiag.6.1.2)), c(sqrt(Lambda_sq.offdiag.6.1.3)), c(sqrt(Lambda_sq.offdiag.6.1.4)), c(sqrt(Lambda_sq.offdiag.6.1.5)), c(sqrt(Lambda_sq.offdiag.6.1.6))), 
                           NuInv=c(c(E_NuInv.offdiag.6.1), c(E_NuInv.offdiag.6.1),  c(E_NuInv.offdiag.6.1),  c(E_NuInv.offdiag.6.1), c(E_NuInv.offdiag.6.1), c(E_NuInv.offdiag.6.1)),
                           truth = factor(c(c(mark.mat), c(mark.mat), c(mark.mat), c(mark.mat), c(mark.mat), c(mark.mat))), 
                           graph = factor(c(rep(1,p^2), rep(2, p^2), rep(3, p^2), rep(4, p^2), rep(5, p^2), rep(6, p^2))))


# Plot theta_ijk's as functions of lambda_ijk's


p1.6.truth <- ggplot2::ggplot(na.omit(df.all.iter.6),  aes(y=theta,x=Lambda))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.6))+
  geom_point(aes(colour=graph, shape=truth))                      

pdf('plots/theta_vs_Lambda_K6_oneeout_100sim_truth.pdf')
p1.6.truth 
dev.off()


# Plot lambda_ijk's as functions of nu_ij's

# Using the same scenario

p2.4.truth <- ggplot2::ggplot(na.omit(df.all.iter.6),  aes(y=Lambda,x=NuInv))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.6))+
  geom_point(aes(colour=graph, shape=truth))                      

pdf('plots/Lambda_vs_NuInv_K6_oneeout_100sim_truth.pdf')
p2.4.truth 
dev.off()

# Plot theta_ijk's as functions of nu_ij's

p3.6.truth <- ggplot2::ggplot(na.omit(df.all.iter.6),  aes(y=theta,x=NuInv))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.6))+
  geom_point(aes(colour=graph, shape=truth))                      

pdf('plots/theta_vs_NuInv_K6_oneeout_100sim_truth.pdf')
p3.6.truth 
dev.off()






