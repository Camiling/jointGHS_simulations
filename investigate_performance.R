source('simulation_functions/perform_jointGHS_simulation.R')
source('simulation_functions/help_functions.R')
library(ggplot2)

# Investigate to behaviour of the jointGHS


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

df.all.iter.1 = data.frame(theta=c(c(theta.offdiag.1.1.1), c(theta.offdiag.1.1.2)), Lambda = c(c(sqrt(Lambda_sq.offdiag.1.1.1)), c(sqrt(Lambda_sq.offdiag.1.1.2))), 
                           NuInv=c(c(E_NuInv.offdiag.1.1), c(E_NuInv.offdiag.1.1)),graph = factor(c(rep(1,p^2), rep(2, p^2))))


# Plot theta_ijk's as functions of lambda_ijk's


p1 <- ggplot2::ggplot(na.omit(df.all.iter.1), aes(y=theta,x=Lambda))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
                   geom_point(aes(colour=graph))                      

pdf('plots/theta_vs_Lambda_K2_80sim.pdf')
p1
dev.off()


# Plot lambda_ijk's as functions of nu_ij's

# Using the same scenario

p2 <- ggplot2::ggplot(na.omit(df.all.iter.1), aes(y=Lambda,x=NuInv))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=graph))                      

pdf('plots/Lambda_vs_NuInv_K2_80sim.pdf')
p2
dev.off()

# Plot theta_ijk's as functions of nu_ij's

p3 <- ggplot2::ggplot(na.omit(df.all.iter.1), aes(y=theta,x=NuInv))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=graph))                      

pdf('plots/theta_vs_NuInv_K2_80sim.pdf')
p3
dev.off()



# Case 2: K=2, identical distributions ----------------------

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

df.all.iter.2 = data.frame(theta=c(c(theta.offdiag.2.1.1), c(theta.offdiag.2.1.2)), Lambda = c(c(sqrt(Lambda_sq.offdiag.2.1.1)), c(sqrt(Lambda_sq.offdiag.2.1.2))), 
                           NuInv=c(c(E_NuInv.offdiag.2.1), c(E_NuInv.offdiag.2.1)),graph = factor(c(rep(1,p^2), rep(2, p^2))))


# Plot theta_ijk's as functions of lambda_ijk's


p1.2 <- ggplot2::ggplot(na.omit(df.all.iter.2), aes(y=theta,x=Lambda))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=graph))                      

pdf('plots/theta_vs_Lambda_K2_100sim.pdf')
p1.2
dev.off()


# Plot lambda_ijk's as functions of nu_ij's

# Using the same scenario

p2.2 <- ggplot2::ggplot(na.omit(df.all.iter.2), aes(y=Lambda,x=NuInv))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=graph))                      

pdf('plots/Lambda_vs_NuInv_K2_100sim.pdf')
p2.2
dev.off()

# Plot theta_ijk's as functions of nu_ij's

p3.2 <- ggplot2::ggplot(na.omit(df.all.iter.2), aes(y=theta,x=NuInv))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour=graph))                      

pdf('plots/theta_vs_NuInv_K2_100sim.pdf')
p3.2
dev.off()






