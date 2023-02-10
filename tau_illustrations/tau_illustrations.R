source('simulation_functions/perform_jointGHS_simulation.R')
source('simulation_functions/help_functions.R')
library(ggplot2)
library(network)
library(GGally)
library(JoStARS)

# Make illustrations of the behavior of tau

# Create trace plots ------------------------------------------

# First example: p=100, n=150

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
tau_sq_vals = seq(0.0001,20,length.out = 500)
theta.est.all = array(0,c(p,p,length(tau_sq_vals)))
for(i in 1:length(tau_sq_vals)){
  res <- fastGHS::fastGHS(x.sf.scaled,tau_sq = tau_sq_vals[i],epsilon = 1e-3, fix_tau=TRUE,verbose = F )
  theta.est = res$theta
  theta.est[which(abs(theta.est)<1e-5, arr.ind=T)]=0
  prec.mat = cov2cor(theta.est)
  theta.est.all[,,i] = prec.mat
  cat(i, '\n')
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

plot.trace <- ggplot2::ggplot(na.omit(trace.df),  aes(y=theta,x=tau_sq, group=edge))+ labs(title=" ")+theme_bw()+theme(plot.title = element_text(hjust = 0.5),text = element_text(size = 15))+
  geom_line(aes(colour=truth, group=edge, linetype=truth))+scale_linetype_manual(values=c("dashed", "solid"))                  

pdf('plots/traceplot_fastGHS_1.pdf')
plot.trace
dev.off()

# Second example: p=150, n=200

n.2 = 200
p.2 = 150
set.seed(123)
data.sf.2 = huge::huge.generator(n=n.2, d=p.2,graph = 'scale-free') 
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
tau_sq_vals.2 = seq(0.0001,20,length.out = 500)
theta.est.all.2 = array(0,c(p.2,p.2,length(tau_sq_vals.2)))
for(i in 1:length(tau_sq_vals.2)){
  res <- fastGHS::fastGHS(x.sf.scaled.2,tau_sq = tau_sq_vals.2[i],epsilon = 1e-3, fix_tau=TRUE, verbose = F)
  theta.est = res$theta
  theta.est[which(abs(theta.est)<1e-5, arr.ind=T)]=0
  prec.mat = cov2cor(theta.est)
  theta.est.all.2[,,i] = prec.mat
  cat(i, '\n')
}

# Create trace plot
theta.est.upper.2 = theta.est.all.2
for(i in 1:length(tau_sq_vals.2)){
  mat.temp = theta.est.all.2[,,i]
  mat.temp[!upper.tri(theta.est.upper.2[,,i])] = NA
  theta.est.upper.2[,,i] = mat.temp
}

# Add information on whether an edge is present in the true network or not
mark.mat.2 = matrix('false',p.2, p.2)
which.mat.2 = theta.true.2!=0
mark.mat.2[which.mat.2==1] = 'true'

trace.df.2 = data.frame(theta=c(theta.est.upper.2), tau_sq=sort(rep(tau_sq_vals.2, length(mark.mat.2))), truth=factor(rep(c(mark.mat.2), length(tau_sq_vals.2))), 
                        edge= factor(rep(1:length(mark.mat.2), length(tau_sq_vals.2))))

plot.trace.2 <- ggplot2::ggplot(na.omit(trace.df.2),  aes(y=theta,x=tau_sq, group=edge))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_line(aes(colour=truth, group=edge, linetype=truth))+scale_linetype_manual(values=c("dashed", "solid"))                          

pdf('plots/traceplot_fastGHS_2.pdf')
plot.trace.2
dev.off()

# Create combined plot

# Get common legend
tmp = ggplot2::ggplot_gtable(ggplot2::ggplot_build(plot.trace+ theme(legend.position = "right")))
leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
legend = tmp$grobs[[leg]]
# Remove legends
plot.trace = plot.trace + theme(legend.position = "none") + labs(title='(a)')
plot.trace.2 = plot.trace.2 + theme(legend.position = "none")+ labs(title='(b)')

pdf(file='plots/traceplot_fastGHS_combined.pdf', 15,8)
gridExtra::grid.arrange(plot.trace,plot.trace.2,legend,ncol=3,widths=c(1,1,0.1))
dev.off()


# Third example: larger n

n.3 = 500
p.3 = 150
set.seed(123)
data.sf.3 = huge::huge.generator(n=n.3, d=p.3,graph = 'scale-free') 
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
tau_sq_vals.3 = seq(0.0001,10,length.out = 200)
theta.est.all.3 = array(0,c(p.3,p.3,length(tau_sq_vals.3)))
for(i in 1:length(tau_sq_vals.3)){
  res <- fastGHS::fastGHS(x.sf.scaled.3,tau_sq = tau_sq_vals.3[i],epsilon = 1e-3, fix_tau=TRUE, verbose = F)
  theta.est = res$theta
  theta.est[which(abs(theta.est)<1e-5, arr.ind=T)]=0
  prec.mat = cov2cor(theta.est)
  theta.est.all.3[,,i] = prec.mat
  cat(i, '\n')
}

# Create trace plot
theta.est.upper.3 = theta.est.all.3
for(i in 1:length(tau_sq_vals.3)){
  mat.temp = theta.est.all.3[,,i]
  mat.temp[!upper.tri(theta.est.upper.3[,,i])] = NA
  theta.est.upper.3[,,i] = mat.temp
}

# Add information on whether an edge is present in the true network or not
mark.mat.3 = matrix('false',p.3, p.3)
which.mat.3 = theta.true.3!=0
mark.mat.3[which.mat.3==1] = 'true'

trace.df.3 = data.frame(theta=c(theta.est.upper.3), tau_sq=sort(rep(tau_sq_vals.3, length(mark.mat.3))), truth=factor(rep(c(mark.mat.3), length(tau_sq_vals.3))), 
                        edge= factor(rep(1:length(mark.mat.3), length(tau_sq_vals.3))))

plot.trace.3 <- ggplot2::ggplot(na.omit(trace.df.3),  aes(y=theta,x=tau_sq, group=edge))+ labs(title=" ")+theme_bw()+theme(plot.title = element_text(hjust = 0.5),text = element_text(size = 15))+
  geom_line(aes(colour=truth, group=edge, linetype=truth))+scale_linetype_manual(values=c("dashed", "solid"))                          

pdf('plots/traceplot_fastGHS_3.pdf')
plot.trace.3
dev.off()

plot.trace.3.2 <- ggplot2::ggplot(na.omit(trace.df.3),  aes(y=theta,x=tau_sq, group=edge))+ labs(title=" ")+theme_bw()+theme(plot.title = element_text(hjust = 0.5),text = element_text(size = 15))+
  geom_line(aes(colour=truth, group=edge, linetype=truth))+scale_linetype_manual(values=c("dashed", "solid")) +scale_x_continuous(trans='log10',breaks = c(0,0.001, 0.01,0.1, 1))                         

pdf('plots/traceplot_fastGHS_3_logscale.pdf')
plot.trace.3.2 
dev.off()

# More zoomed in version:

# Perform GHS for a range of tau_sq values
tau_sq_vals.3.2 = c(seq(0.0001,0.015,length.out = 200),seq(0.015,0.15,length.out = 200))
theta.est.all.3.2 = array(0,c(p.3,p.3,length(tau_sq_vals.3.2)))
for(i in 1:length(tau_sq_vals.3.2)){
  res <- fastGHS::fastGHS(x.sf.scaled.3,tau_sq = tau_sq_vals.3.2[i],epsilon = 1e-3, fix_tau=TRUE, verbose = F)
  theta.est = res$theta
  theta.est[which(abs(theta.est)<1e-5, arr.ind=T)]=0
  prec.mat = cov2cor(theta.est)
  theta.est.all.3.2[,,i] = prec.mat
  cat(i, '\n')
}

# Create trace plot
theta.est.upper.3.2 = theta.est.all.3.2
for(i in 1:length(tau_sq_vals.3.2)){
  mat.temp = theta.est.all.3.2[,,i]
  mat.temp[!upper.tri(theta.est.upper.3.2[,,i])] = NA
  theta.est.upper.3.2[,,i] = mat.temp
}

trace.df.3.2 = data.frame(theta=c(theta.est.upper.3.2), tau_sq=sort(rep(tau_sq_vals.3.2, length(mark.mat.3))), truth=factor(rep(c(mark.mat.3), length(tau_sq_vals.3.2))), 
                        edge= factor(rep(1:length(mark.mat.3), length(tau_sq_vals.3.2))))

plot.trace.3.2.zoom <- ggplot2::ggplot(na.omit(trace.df.3.2),  aes(y=theta,x=tau_sq, group=edge))+ labs(title=" ")+theme_bw()+theme(plot.title = element_text(hjust = 0.5),text = element_text(size = 15))+
  geom_line(aes(colour=truth, group=edge, linetype=truth))+scale_linetype_manual(values=c("dashed", "solid")) +scale_x_continuous(trans='log10',breaks = c(0,0.001, 0.01,0.1, 1))                         

pdf('plots/traceplot_fastGHS_3_logscale_smaller.pdf')
plot.trace.3.2.zoom
dev.off()

# Combined plot
library(patchwork)

pdf('plots/traceplot_comb.pdf', 10,5)
plot.trace.3.2.zoom + theme(legend.position="none")+ ggtitle('(a)') + (plot.trace.3 + ggtitle('(b)')) 
dev.off()





# Plot tau vs sparsity, precision and recall  ---------------------------------------------------------

# New example
p=100
n=100
set.seed(12)
data.sf = huge::huge.generator(n=n, d=p,graph = 'scale-free') 
g.true.sf = data.sf$theta # True adjacency matrix
theta.true = data.sf$omega # The precision matrix
theta.true[which(theta.true<10e-5,arr.ind=T)]=0  
g.sf=graph.adjacency(data.sf$theta,mode="undirected",diag=F) # true igraph object
x.sf = data.sf$data # Observed attributes. nxp matrix.
x.sf.scaled= scale(x.sf) # Scale columns/variables.s.sf.scaled = cov(x.sf.scaled) # Empirical covariance matrix
s.sf.scaled = cov(x.sf.scaled) # Empirical covariance matrix


# Perform GHS for a range of tau_sq values, first network
tau_sq_vals.1.1 = seq(1e-2, 10, length.out=50)
theta.est.all.1.1 = array(0,c(p,p,length(tau_sq_vals.1.1)))
for(i in 1:length(tau_sq_vals.1.1)){
  res <- fastGHS::fastGHS(x.sf.scaled,tau_sq = tau_sq_vals.1.1[i],epsilon = 1e-3, fix_tau=TRUE,verbose = F)
  theta.est = res$theta
  theta.est[which(abs(theta.est)<1e-5, arr.ind=T)]=0
  prec.mat = cov2cor(theta.est)
  theta.est.all.1.1[,,i] = prec.mat
  cat(i, '\n')
}

tau_sq_vals.1.2 = seq(1e-2, 10, length.out=50)
theta.est.all.1.2 = array(0,c(p.2,p.2,length(tau_sq_vals.1.2)))
for(i in 1:length(tau_sq_vals.1.2)){
  res <- fastGHS::fastGHS(x.sf.scaled.2,tau_sq = tau_sq_vals.1.2[i],epsilon = 1e-3, fix_tau=TRUE,verbose = F)
  theta.est = res$theta
  theta.est[which(abs(theta.est)<1e-5, arr.ind=T)]=0
  prec.mat = cov2cor(theta.est)
  theta.est.all.1.2[,,i] = prec.mat
  cat(i, '\n')
}

# Sparsities
spars.1= unlist(apply(theta.est.all.1.1, 3, FUN = function(s) tailoredGlasso::sparsity(s!=0)))
spars.2 = unlist(apply(theta.est.all.1.2, 3, FUN = function(s) tailoredGlasso::sparsity(s!=0)))

# Precisions 
prec.1 = unlist(apply(theta.est.all.1.1, 3, FUN = function(s) tailoredGlasso::precision(theta.true!=0, s!=0)))
prec.2 = unlist(apply(theta.est.all.1.2, 3,  FUN = function(s) tailoredGlasso::precision(theta.true.2!=0, s!=0)))

# Recalls 
rec.1 = unlist(apply(theta.est.all.1.1, 3, FUN = function(s) tailoredGlasso::recall(theta.true!=0, s!=0)))
rec.2 = unlist(apply(theta.est.all.1.2, 3, FUN = function(s) tailoredGlasso::recall(theta.true.2!=0, s!=0)))

df.1 = data.frame(tau_sq = tau_sq_vals.1.1, sparsity = spars.1)
df.2 = data.frame(tau_sq = tau_sq_vals.1.2, sparsity = spars.2)
df.prec.1 = data.frame(tau_sq = tau_sq_vals.1.1, precision = prec.1)
df.prec.2 = data.frame(tau_sq = tau_sq_vals.1.2, precision = prec.2)
df.rec.1 = data.frame(tau_sq = tau_sq_vals.1.1, recall = rec.1)
df.rec.2 = data.frame(tau_sq = tau_sq_vals.1.2, recall = rec.2)


p1 <- ggplot2::ggplot(df.1, aes(x=tau_sq,y=sparsity))+ geom_line(color='steelblue', size=1)+ labs(title="(a)")+theme_bw()+theme(plot.title = element_text(hjust = 0.5),text = element_text(size = 15))+
          scale_y_continuous(limits=c(0,0.015))
p2 <- ggplot2::ggplot(df.2, aes(x=tau_sq,y=sparsity))+ geom_line(color='steelblue', size=1)+ labs(title="(b)")+theme_bw()+theme(plot.title = element_text(hjust = 0.5),text = element_text(size = 15))+
          scale_y_continuous(limits=c(0,0.015))
p1.prec <- ggplot2::ggplot(df.prec.1, aes(x=tau_sq,y=precision))+ geom_line(color='steelblue', size=1)+ labs(title=" ")+theme_bw()+theme(plot.title = element_text(hjust = 0.5),text = element_text(size = 15))+
            scale_y_continuous(limits=c(0,1))
p2.prec <- ggplot2::ggplot(df.prec.2, aes(x=tau_sq,y=precision))+ geom_line(color='steelblue', size=1)+ labs(title=" ")+theme_bw()+theme(plot.title = element_text(hjust = 0.5),text = element_text(size = 15))+
            scale_y_continuous(limits=c(0,1))
p1.rec <- ggplot2::ggplot(df.rec.1, aes(x=tau_sq,y=recall))+ geom_line(color='steelblue', size=1)+ labs(title=" ")+theme_bw()+theme(plot.title = element_text(hjust = 0.5),text = element_text(size = 15))+
            scale_y_continuous(limits=c(0,1))
p2.rec <- ggplot2::ggplot(df.rec.2, aes(x=tau_sq,y=recall))+ geom_line(color='steelblue', size=1)+ labs(title=" ")+theme_bw()+theme(plot.title = element_text(hjust = 0.5),text = element_text(size = 15))+
            scale_y_continuous(limits=c(0,1))

# Plot sparsity, precision and recall as functions of tau

pdf('plots/measures_vs_tau.pdf')
gridExtra::grid.arrange(p1, p2, p1.prec, p2.prec, p1.rec, p2.rec,ncol=2)
dev.off()


# Plot tau vs AIC  --------------------------------------------------------------

# Computee the AIC for different tau_sq

tau_sq_vals.aic = seq(0.001,2,length.out = 100)
sparsities.aic = rep(0,length(tau_sq_vals.aic))
prec.aic = rep(0,length(tau_sq_vals.aic))
rec.aic = rep(0,length(tau_sq_vals.aic))
aic.vals = rep(0,length(tau_sq_vals.aic))
for(i in 1:length(tau_sq_vals.aic)){
  res <- fastGHS::fastGHS(x.sf.scaled,tau_sq = tau_sq_vals.aic[i],epsilon = 1e-3, fix_tau=TRUE, verbose=F)
  theta.est = res$theta
  theta.est[which(abs(theta.est)<1e-5, arr.ind=T)]=0
  prec.mat = cov2cor(theta.est)
  aic.vals[i] = gaussianAIC(s.sf.scaled,prec.mat, n)
  sparsities.aic[i] = tailoredGlasso::sparsity(theta.est!=0)
  prec.aic[i] = tailoredGlasso::precision(theta.true!=0, theta.est!=0)
  rec.aic[i] = tailoredGlasso::recall(theta.true!=0, theta.est!=0)
  cat(i,'\n')
}

tau_sq_vals.aic.2 = seq(0.001,2,length.out = 100)
sparsities.aic.2 = rep(0,length(tau_sq_vals.aic.2))
prec.aic.2 = rep(0,length(tau_sq_vals.aic.2))
rec.aic.2 = rep(0,length(tau_sq_vals.aic.2))
aic.vals.2 = rep(0,length(tau_sq_vals.aic.2))
for(i in 1:length(tau_sq_vals.aic.2)){
  res <- fastGHS::fastGHS(x.sf.scaled.2,tau_sq = tau_sq_vals.aic.2[i],epsilon = 1e-3, fix_tau=TRUE, verbose=F)
  theta.est = res$theta
  theta.est[which(abs(theta.est)<1e-5, arr.ind=T)]=0
  prec.mat = cov2cor(theta.est)
  aic.vals.2[i] = gaussianAIC(s.sf.scaled.2,prec.mat, n.2)
  sparsities.aic.2[i] = tailoredGlasso::sparsity(theta.est!=0)
  prec.aic.2[i] = tailoredGlasso::precision(theta.true.2!=0, theta.est!=0)
  rec.aic.2[i] = tailoredGlasso::recall(theta.true.2!=0, theta.est!=0)
  cat(i,'\n')
}

df.aic = data.frame(sparsity=sparsities.aic[-1], AIC = aic.vals[-1], precision = prec.aic[-1], recall=rec.aic[-1], tau_sq=tau_sq_vals.aic[-1])
df.aic.2 = data.frame(sparsity=sparsities.aic.2[-1], AIC = aic.vals.2[-1], precision = prec.aic.2[-1], recall=rec.aic.2[-1], tau_sq=tau_sq_vals.aic.2[-1])
which.conv = min(which(abs(df.aic$AIC[2:length(df.aic$AIC)]-df.aic$AIC[1:(length(df.aic$AIC)-1)]) < 0.1))+1 # Add one to get first value with AIC diff less than 0.1 from the previous
which.conv.2 = min(which(abs(df.aic.2$AIC[2:length(df.aic.2$AIC)]-df.aic.2$AIC[1:(length(df.aic.2$AIC)-1)]) < 0.1))+1
tau_sq_opt = tau_sq_vals.aic[which.conv]
tau_sq_opt.2 = tau_sq_vals.aic.2[which.conv.2]


plot.aic = ggplot2::ggplot(df.aic,  aes(y=AIC,x=tau_sq))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_line(color='steelblue', size=1) +geom_point(aes(x=tau_sq_opt,aic.vals[which.conv]), color='hotpink1', size=2)+ labs(title="(a)")
plot.spars = ggplot2::ggplot(df.aic,  aes(y=sparsity,x=tau_sq))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_line(colour='steelblue', size=1)  +geom_point(aes(x=tau_sq_opt,sparsities.aic[which.conv]), color='hotpink1', size=2)  +scale_y_continuous(limits=c(0,0.015))
plot.prec = ggplot2::ggplot(df.aic,  aes(y=precision,x=tau_sq))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_line(colour='steelblue', size=1)  +geom_point(aes(x=tau_sq_opt,prec.aic[which.conv]), color='hotpink1', size=2)  +scale_y_continuous(limits=c(0,1))
plot.rec = ggplot2::ggplot(df.aic,  aes(y=recall,x=tau_sq))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_line(colour='steelblue', size=1) +geom_point(aes(x=tau_sq_opt,rec.aic[which.conv]), color='hotpink1', size=2)  +scale_y_continuous(limits=c(0,1))

plot.aic.2 = ggplot2::ggplot(df.aic.2,  aes(y=AIC,x=tau_sq))+ labs(title=" ")+theme_bw()+theme(plot.title = element_text(hjust = 0.5),text = element_text(size = 15))+
  geom_line(color='steelblue', size=1) +geom_point(aes(x=tau_sq_opt.2,aic.vals.2[which.conv.2]), color='hotpink1', size=2)+ labs(title="(b)")
plot.spars.2 = ggplot2::ggplot(df.aic.2,  aes(y=sparsity,x=tau_sq))+ labs(title=" ")+theme_bw()+theme(plot.title = element_text(hjust = 0.5),text = element_text(size = 15))+
  geom_line(colour='steelblue', size=1) +geom_point(aes(x=tau_sq_opt.2,sparsities.aic.2[which.conv.2]), color='hotpink1', size=2)+scale_y_continuous(limits=c(0,0.015))
plot.prec.2 = ggplot2::ggplot(df.aic.2,  aes(y=precision,x=tau_sq))+ labs(title=" ")+theme_bw()+theme(plot.title = element_text(hjust = 0.5),text = element_text(size = 15))+
  geom_line(colour='steelblue', size=1)  +geom_point(aes(x=tau_sq_opt.2,prec.aic.2[which.conv.2]), color='hotpink1', size=2) +scale_y_continuous(limits=c(0,1))
plot.rec.2 = ggplot2::ggplot(df.aic.2,  aes(y=recall,x=tau_sq))+ labs(title=" ")+theme_bw()+theme(plot.title = element_text(hjust = 0.5),text = element_text(size = 15))+
  geom_line(colour='steelblue', size=1) +geom_point(aes(x=tau_sq_opt.2,rec.aic.2[which.conv.2]), color='hotpink1', size=2) +scale_y_continuous(limits=c(0,1))


pdf('plots/AIC_and_measures_vs_tau.pdf', 10,20)
gridExtra::grid.arrange(plot.aic,plot.aic.2, plot.spars,  plot.spars.2, plot.prec, plot.prec.2, plot.rec, plot.rec.2,ncol=2)
dev.off()

# Save only one network
pdf('plots/AIC_and_measures_vs_tau_onlyone.pdf', 10,10)
gridExtra::grid.arrange(plot.aic.2+labs(title=""), plot.spars.2, plot.prec.2,plot.rec.2,ncol=2)
dev.off()

# Plot tau_sq trajectory in ECM algorithm --------------------------------------

# For a range of p and n

p.vals = c(10, 20, 50, 80, 100, 150,200,250, 300, 350, 400, 450, 500)
n.val = 500
tau_sq_all = list()
p.vals.all= list()
iter.no.all = list()

set.seed(123)
for(i in 1:length(p.vals)){
  data.tmp = huge::huge.generator(n=n.val, d=p.vals[i],graph = 'scale-free') 
  g.tmp = data.tmp$theta # True adjacency matrix
  theta.true.tmp = data.tmp$omega # The precision matrix
  theta.true.tmp[which(theta.true.tmp<10e-5,arr.ind=T)]=0  
  x.tmp = data.tmp$data # Observed attributes. nxp matrix.
  x.tmp.scaled= scale(x.tmp) # Scale columns/variables.
  res.tmp = fastGHS::fastGHS(x.tmp.scaled, epsilon = 1e-3, fix_tau = F,AIC_selection = F, verbose=F, savepath_tau = T)
  tau_sq_all[[i]] = res.tmp$tau_sq_path[-1]
  p.vals.all[[i]] = rep(p.vals[i], length(res.tmp$tau_sq_path)-1)
  iter.no.all[[i]] = 1:(length(res.tmp$tau_sq_path)-1)
  cat(i, '\n')
}

df.tau = data.frame(tau_sq = unlist(tau_sq_all), p=factor(unlist(p.vals.all)), iter=unlist(iter.no.all))

p.shrink = ggplot2::ggplot(df.tau, aes(y=tau_sq, x=iter, group=p))+geom_line(aes(col=p))

pdf('plots/tau_shinking.pdf')
p.shrink
dev.off()


# Random initial values
p.init=50
n.init=500
set.seed(123)
data.tmp = huge::huge.generator(n=n.init, d=p.init,graph = 'scale-free') 
g.tmp = data.tmp$theta # True adjacency matrix
theta.true.tmp = data.tmp$omega # The precision matrix
theta.true.tmp[which(theta.true.tmp<10e-5,arr.ind=T)]=0  
x.tmp = data.tmp$data # Observed attributes. nxp matrix.
x.tmp.scaled= scale(x.tmp) # Scale columns/variables.

tau_sq_init_vals = seq(0.1,2, by=0.1)
tau_sq_all.init = list()
init_vals.init= list()
iter.no.all.init = list()
for(i in 1:length(tau_sq_init_vals)){
  res.tmp = fastGHS::fastGHS(x.tmp.scaled, epsilon = 1e-3, tau_sq=tau_sq_init_vals[i], fix_tau = F,AIC_selection = F, verbose=F, savepath_tau = T)
  tau_sq_all.init[[i]] = res.tmp$tau_sq_path[-1]
  init_vals.init[[i]] = rep(tau_sq_init_vals[i], length(res.tmp$tau_sq_path)-1)
  iter.no.all.init[[i]] = 1:(length(res.tmp$tau_sq_path)-1)
  cat(i, '\n')
}

df.tau.init = data.frame(tau_sq = unlist(tau_sq_all.init), tau_sq_init=factor(unlist(init_vals.init)), iter=unlist(iter.no.all.init))
p.init<- ggplot2::ggplot(df.tau.init, aes(y=tau_sq, x=iter, group=tau_sq_init))+geom_line(aes(col=tau_sq_init))+theme_bw()+theme(text = element_text(size = 15)) 
p.init$labels$colour <- "tau_sq initial value"

pdf('plots/tau_shinking_init.pdf')
p.init
dev.off()

# Make combined plots

pdf('plots/tau_shinking_combined.pdf', 18,8)
gridExtra::grid.arrange(p.shrink+labs(title="(a)")+theme(plot.title = element_text(hjust = 0.5)), p.init+labs(title="(b)")+theme(plot.title = element_text(hjust = 0.5)),ncol=2)
dev.off()






