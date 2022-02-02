source('simulation_functions/perform_jointGHS_simulation.R')
source('simulation_functions/help_functions.R')
library(ggplot2)
library(network)
library(GGally)

# Create trace plots

# EXAMPLE 1 ------------------------------------------------------------------
# GENERATE GRAPH: n=150, p=100

n=150
p=100
set.seed(123)
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

tau_sq_vals = seq(0.0001,0.3,length.out = 100)
theta.est.all = array(0,c(p,p,length(tau_sq_vals)))
for(i in 1:length(tau_sq_vals)){
  res <- fastGHS::fastGHS(x.sf.scaled,tau_sq = tau_sq_vals[i],epsilon = 1e-3, fix_tau=TRUE)
  theta.est = res$theta
  theta.est[which(abs(theta.est)<1e-5, arr.ind=T)]=0
  prec.mat = cov2cor(theta.est)
  theta.est.all[,,i] = prec.mat
}


# Create trace plot

theta.est.upper = theta.est.all
for(i in 1:length(tau_sq_vals)){
  mat.temp = theta.est.all[,,i]
  mat.temp[!upper.tri(theta.est.upper[,,i])] = NA
  theta.est.upper[,,i] = mat.temp
}

# Add information on whether an edge is present in the true network or not
mark.mat = matrix('false',p, p)
which.mat = theta.true!=0
mark.mat[which.mat==1] = 'true'

trace.df = data.frame(theta=c(theta.est.upper), tau_sq=sort(rep(tau_sq_vals, length(mark.mat))), truth=factor(rep(c(mark.mat), length(tau_sq_vals))), 
                      edge= factor(rep(1:length(mark.mat), length(tau_sq_vals))))

plot.trace <- ggplot2::ggplot(na.omit(trace.df),  aes(y=theta,x=tau_sq, group=edge))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_line(aes(colour=truth, group=edge))                      

pdf('plots/plots2/traceplot_fastGHS.pdf')
plot.trace
dev.off()

plot.trace.small <- ggplot2::ggplot(na.omit(trace.df),  aes(y=theta,x=tau_sq, group=edge))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_line(aes(colour=truth, group=edge))+xlim(0,0.075)                      


pdf('plots/plots2/traceplot_fastGHS_small.pdf')
plot.trace.small
dev.off()

# Investigate the one negative value

ind.neg = which(theta.est.upper[,,50]< -1e-5,arr.ind = T)
ind.neg
# row col
# [1,]   2  53

# Look at their sample covariance
samplecov = s.sf.scaled[upper.tri(s.sf.scaled)]

pdf('plots/plots2/neg_cov_fastGHS.pdf')
hist(samplecov, breaks=100)
abline(v=s.sf.scaled[2,53], col='red')
dev.off()

# The sample covariance is very large for these two nodes. 

# Check common neighbours:
which(theta.true[2,]!=0 &theta.true[53,]!=0)
# 1

# How big is the partial correlation they have to this common neighbour?
theta.true[2,1]
theta.true[53,1]
partial.cors = cov2cor(theta.true)[upper.tri(theta.true)]
pdf('plots/plots2/negcov_common_fastGHS.pdf')
hist(partial.cors, breaks=200)
abline(v=cov2cor(theta.true)[2,1], col='red')
abline(v=cov2cor(theta.true)[53,1], col='blue')
dev.off()


# Plot graph with these two highlighted
net = network::network(theta.true!=0)
col.net =rep("deepskyblue2",ncol(theta.true))
col.net[c(2,53)] = "hotpink1"

pdf('plots/plots2/graph_fastGHS_negval.pdf')
set.seed(123)
GGally::ggnet2(net, alpha = 0.9, mode = "fruchtermanreingold", color = col.net) 
dev.off()

# Check their lambda value:
res.test = fastGHS(x.sf.scaled, tau_sq=0.1)
res.test$Lambda_sq[2,53]
# 0.04604704

lambdas.upper = res.test$Lambda_sq[upper.tri(res.test$Lambda_sq)]

pdf('plots/plots2/Lambda_fastGHS_negval.pdf')
hist(lambdas.upper, breaks=200)
abline(v=res.test$Lambda_sq[2,53], col='red')
dev.off()


# EXAMPLE 2 ------------------------------------------------------------------
# GENERATE GRAPH: n=200, p=150

n=200
p=150
set.seed(123)
data.sf.2 = huge::huge.generator(n=n, d=p,graph = 'scale-free') 
g.true.sf.2 = data.sf.2$theta # True adjacency matrix
theta.true.2 = data.sf.2$omega # The precision matrix
theta.true.2[which(theta.true.2<10e-5,arr.ind=T)]=0  
g.sf.2=graph.adjacency(data.sf.2$theta,mode="undirected",diag=F) # true igraph object
x.sf.2 = data.sf.2$data # Observed attributes. nxp matrix.
x.sf.scaled.2= scale(x.sf.2) # Scale columns/variables.
s.sf.scaled.2 = cov(x.sf.scaled.2) # Empirical covariance matrix
data.sf.2$sparsity
# 0.01333333


# Perform GHS for a range of tau_sq values

tau_sq_vals.2 = seq(0.0001,2,length.out = 100)
theta.est.all.2 = array(0,c(p,p,length(tau_sq_vals.2)))
for(i in 1:length(tau_sq_vals.2)){
  res <- fastGHS::fastGHS(x.sf.scaled.2,tau_sq = tau_sq_vals.2[i],epsilon = 1e-3, fix_tau=TRUE)
  theta.est = res$theta
  theta.est[which(abs(theta.est)<1e-5, arr.ind=T)]=0
  prec.mat = cov2cor(theta.est)
  theta.est.all.2[,,i] = prec.mat
}


# Create trace plot

theta.est.upper.2 = theta.est.all.2
for(i in 1:length(tau_sq_vals.2)){
  mat.temp = theta.est.all.2[,,i]
  mat.temp[!upper.tri(theta.est.upper.2[,,i])] = NA
  theta.est.upper.2[,,i] = mat.temp
}

# Add information on whether an edge is present in the true network or not
mark.mat = matrix('false',p, p)
which.mat = theta.true.2!=0
mark.mat[which.mat==1] = 'true'

trace.df.2 = data.frame(theta=c(theta.est.upper.2), tau_sq=sort(rep(tau_sq_vals.2, length(mark.mat))), truth=factor(rep(c(mark.mat), length(tau_sq_vals.2))), 
                      edge= factor(rep(1:length(mark.mat), length(tau_sq_vals.2))))

plot.trace.2 <- ggplot2::ggplot(na.omit(trace.df.2),  aes(y=theta,x=tau_sq, group=edge))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_line(aes(colour=truth, group=edge))                      

pdf('plots/plots2/traceplot_fastGHS_2.pdf')
plot.trace.2
dev.off()

plot.trace.small.2 <- ggplot2::ggplot(na.omit(trace.df.2),  aes(y=theta,x=tau_sq, group=edge))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_line(aes(colour=truth, group=edge))+xlim(0,0.075)                      


pdf('plots/plots2/traceplot_fastGHS_small_2.pdf')
plot.trace.small.2
dev.off()

# EXAMPLE 3: Graph 2, but with larger n ------------------------------------------------------------------
# GENERATE GRAPH: n=500, p=150

n=500
p=150
set.seed(123)
data.sf.3 = huge::huge.generator(n=n, d=p,graph = 'scale-free') 
g.true.sf.3 = data.sf.3$theta # True adjacency matrix
theta.true.3 = data.sf.3$omega # The precision matrix
theta.true.3[which(theta.true.3<10e-5,arr.ind=T)]=0  
g.sf.3=graph.adjacency(data.sf.3$theta,mode="undirected",diag=F) # true igraph object
x.sf.3 = data.sf.3$data # Observed attributes. nxp matrix.
x.sf.scaled.3= scale(x.sf.3) # Scale columns/variables.
s.sf.scaled.3 = cov(x.sf.scaled.3) # Empirical covariance matrix
data.sf.3$sparsity
# 0.02


# Perform GHS for a range of tau_sq values

tau_sq_vals.3 = seq(0.0001,2,length.out = 100)
theta.est.all.3 = array(0,c(p,p,length(tau_sq_vals.3)))
for(i in 1:length(tau_sq_vals.3)){
  res <- fastGHS::fastGHS(x.sf.scaled.3,tau_sq = tau_sq_vals.3[i],epsilon = 1e-3, fix_tau=TRUE)
  theta.est = res$theta
  theta.est[which(abs(theta.est)<1e-5, arr.ind=T)]=0
  prec.mat = cov2cor(theta.est)
  theta.est.all.3[,,i] = prec.mat
}


# Create trace plot

theta.est.upper.3 = theta.est.all.3
for(i in 1:length(tau_sq_vals.3)){
  mat.temp = theta.est.all.3[,,i]
  mat.temp[!upper.tri(theta.est.upper.3[,,i])] = NA
  theta.est.upper.3[,,i] = mat.temp
}

# Add information on whether an edge is present in the true network or not
mark.mat = matrix('false',p, p)
which.mat = theta.true.3!=0
mark.mat[which.mat==1] = 'true'

trace.df.3 = data.frame(theta=c(theta.est.upper.3), tau_sq=sort(rep(tau_sq_vals.3, length(mark.mat))), truth=factor(rep(c(mark.mat), length(tau_sq_vals.3))), 
                        edge= factor(rep(1:length(mark.mat), length(tau_sq_vals.3))))

plot.trace.3 <- ggplot2::ggplot(na.omit(trace.df.3),  aes(y=theta,x=tau_sq, group=edge))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_line(aes(colour=truth, group=edge))                      

pdf('plots/plots2/traceplot_fastGHS_3.pdf')
plot.trace.3
dev.off()

plot.trace.small.3 <- ggplot2::ggplot(na.omit(trace.df.3),  aes(y=theta,x=tau_sq, group=edge))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_line(aes(colour=truth, group=edge))+xlim(0,0.075)                      


pdf('plots/plots2/traceplot_fastGHS_small_3.pdf')
plot.trace.small.3
dev.off()

# EXAMPLE 4: Graph 1, but with a finer grid around 0 ----------------------------------

n=150
p=100
set.seed(123)
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

tau_sq_vals = c(seq(0.0001,0.02,length.out = 300), seq(0.0201,0.3,length.out = 100))
theta.est.all = array(0,c(p,p,length(tau_sq_vals)))
for(i in 1:length(tau_sq_vals)){
  res <- fastGHS::fastGHS(x.sf.scaled,tau_sq = tau_sq_vals[i],epsilon = 1e-3, fix_tau=TRUE)
  theta.est = res$theta
  theta.est[which(abs(theta.est)<1e-5, arr.ind=T)]=0
  prec.mat = cov2cor(theta.est)
  theta.est.all[,,i] = prec.mat
}


# Create trace plot

theta.est.upper = theta.est.all
for(i in 1:length(tau_sq_vals)){
  mat.temp = theta.est.all[,,i]
  mat.temp[!upper.tri(theta.est.upper[,,i])] = NA
  theta.est.upper[,,i] = mat.temp
}

# Add information on whether an edge is present in the true network or not
mark.mat = matrix('false',p, p)
which.mat = theta.true!=0
mark.mat[which.mat==1] = 'true'

trace.df = data.frame(theta=c(theta.est.upper), tau_sq=sort(rep(tau_sq_vals, length(mark.mat))), truth=factor(rep(c(mark.mat), length(tau_sq_vals))), 
                      edge= factor(rep(1:length(mark.mat), length(tau_sq_vals))))

plot.trace <- ggplot2::ggplot(na.omit(trace.df),  aes(y=theta,x=tau_sq, group=edge))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_line(aes(colour=truth, group=edge))                      

pdf('plots/plots2/traceplot_fastGHS_finergrid.pdf')
plot.trace
dev.off()

plot.trace.small <- ggplot2::ggplot(na.omit(trace.df),  aes(y=theta,x=tau_sq, group=edge))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_line(aes(colour=truth, group=edge))+xlim(0,0.075)                      


pdf('plots/plots2/traceplot_fastGHS_finergrid_small.pdf')
plot.trace.small
dev.off()


