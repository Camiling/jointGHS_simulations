rm(list=ls())
source('simulation_functions/perform_jointGHS_simulation.R')
source('simulation_functions/help_functions.R')


# Plot results from single vs joint GHS

load("data/jointGHS_simulations_K2_singleVSjoint_final.Rdata")
load("data/jointGHS_simulations_K4_singleVSjoint_final.Rdata") 
load("data/jointGHS_simulations_K10_singleVSjoint_final.Rdata") 


fracs.disagreement = c(0,0.2,0.4,0.6,0.8,1)
perc.disagreement = 100*fracs.disagreement

precisions.K2 = unlist(lapply(res.K2, FUN = function(s) mean(s$mean.precisions[1])))
precisions.K2.ghs = unlist(lapply(res.K2, FUN = function(s) mean(s$mean.precisions.ghs[1])))

precisions.K4 = unlist(lapply(res.K4.list, FUN = function(s) mean(s$mean.precisions[1])))
precisions.K4.ghs = unlist(lapply(res.K4.list, FUN = function(s) mean(s$mean.precisions.ghs[1])))

precisions.K10 = unlist(lapply(res.K10.list, FUN = function(s) mean(s$mean.precisions[1])))
precisions.K10.ghs = unlist(lapply(res.K10.list, FUN = function(s) mean(s$mean.precisions.ghs[1])))

recalls.K2 = unlist(lapply(res.K2, FUN = function(s) mean(s$mean.recalls[1])))
recalls.K2.ghs = unlist(lapply(res.K2, FUN = function(s) mean(s$mean.recalls.ghs[1])))

recalls.K4 = unlist(lapply(res.K4.list, FUN = function(s) mean(s$mean.recalls[1])))
recalls.K4.ghs = unlist(lapply(res.K4.list, FUN = function(s) mean(s$mean.recalls.ghs[1])))

recalls.K10 = unlist(lapply(res.K10.list, FUN = function(s) mean(s$mean.recalls[1])))
recalls.K10.ghs = unlist(lapply(res.K10.list, FUN = function(s) mean(s$mean.recalls.ghs[1])))

spars.K2 = unlist(lapply(res.K2, FUN = function(s) mean(s$mean.opt.sparsities[1])))
spars.K2.ghs = unlist(lapply(res.K2, FUN = function(s) mean(s$mean.opt.sparsities.ghs[1])))

spars.K4= unlist(lapply(res.K4.list, FUN = function(s) mean(s$mean.opt.sparsities[1])))
spars.K4.ghs = unlist(lapply(res.K4.list, FUN = function(s) mean(s$mean.opt.sparsities.ghs[1])))

spars.K10= unlist(lapply(res.K10.list, FUN = function(s) mean(s$mean.opt.sparsities[1])))
spars.K10.ghs = unlist(lapply(res.K10.list, FUN = function(s) mean(s$mean.opt.sparsities.ghs[1])))

df.K = data.frame(precision=c(precisions.K2, precisions.K2.ghs,precisions.K4,precisions.K4.ghs,precisions.K10,precisions.K10.ghs), 
                  recall=c(recalls.K2, recalls.K2.ghs,recalls.K4,recalls.K4.ghs,recalls.K10,recalls.K10.ghs),
                  sparsity=c(spars.K2, spars.K2.ghs,spars.K4,spars.K4.ghs,spars.K10,spars.K10.ghs),
                  method = factor(c(rep('jointGHS', length(precisions.K2)),rep('fastGHS',length(precisions.K2.ghs)),
                                         rep('jointGHS', length(precisions.K4)),rep('fastGHS',length(precisions.K4.ghs)),
                                         rep('jointGHS', length(precisions.K10)),rep('fastGHS',length(precisions.K10.ghs))), levels=c('jointGHS','fastGHS')),
                       K = factor(c(rep('K=2', length(precisions.K2)),rep('K=2',length(precisions.K2.ghs)),
                                    rep('K=4', length(precisions.K4)),rep('K=4',length(precisions.K4.ghs)),
                                    rep('K=10', length(precisions.K10)),rep('K=10',length(precisions.K10.ghs))),levels=c('K=2','K=4','K=10')),
                       disagreement=c(rep(perc.disagreement,6)))

g.prec=ggplot2::ggplot(df.K, aes(y=precision, x=disagreement, group=method))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_line(aes(colour=method, linetype=method), size=1)+ scale_color_manual(values=c("darkolivegreen", "darkolivegreen3"))+
  labs(x='Disagreement %')+ scale_x_continuous(breaks = perc.disagreement) + facet_wrap(~K) + theme(legend.position = "bottom")#+ scale_y_continuous(breaks = seq(0,60,by=10))


g.rec = ggplot2::ggplot(df.K, aes(y=recall, x=disagreement, group=method))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_line(aes(colour=method, linetype=method), size=1)+ scale_color_manual(values=c("darkolivegreen", "darkolivegreen3"))+
  labs(x='Disagreement %')+ scale_x_continuous(breaks = perc.disagreement) + facet_wrap(~K)#+ scale_y_continuous(breaks = seq(0,60,by=10))

g.spars= ggplot2::ggplot(df.K, aes(y=sparsity, x=disagreement, group=method))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_line(aes(colour=method, linetype=method), size=1)+ scale_color_manual(values=c("darkolivegreen", "darkolivegreen3"))+
  labs(x='Disagreement %')+ scale_x_continuous(breaks = perc.disagreement) + facet_wrap(~K)#+ scale_y_continuous(breaks = seq(0,60,by=10))


tmp = ggplot2::ggplot_gtable(ggplot2::ggplot_build(g.prec))
leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
legend = tmp$grobs[[leg]]
g.rec = g.rec + theme_bw()+theme(legend.position = "none",text = element_text(size = 15))
g.prec = g.prec + theme_bw()+theme(legend.position = "none",text = element_text(size = 15))
g.spars = g.spars+ theme_bw()+theme(legend.position = "none",text = element_text(size = 15))

pdf(file='plots/joint_vs_single.pdf',10,8)
gridExtra::grid.arrange(g.prec,g.rec,legend,ncol=1,heights=c(6,6,1))
dev.off()

pdf(file='plots/joint_vs_single_withspars.pdf',10,10)
gridExtra::grid.arrange(g.prec,g.rec,g.spars,legend,ncol=1,heights=c(6,6,6,1))
dev.off()

# True sparsity 
res.K2[[1]]$true.sparsity



# Make plots with SE shown as well

# Only for the first network
precisions.sd.K2 = unlist(lapply(res.K2, FUN = function(s) sd(s$precisions[,1])))
precisions.ghs.sd.K2 = unlist(lapply(res.K2, FUN = function(s) sd(s$precisions.ghs[,1])))

precisions.sd.K4 = unlist(lapply(res.K4.list, FUN = function(s) sd(s$precisions[,1])))
precisions.ghs.sd.K4 = unlist(lapply(res.K4.list, FUN = function(s) sd(s$precisions.ghs[,1])))

precisions.sd.K10 = unlist(lapply(res.K10.list, FUN = function(s) sd(s$precisions[,1])))
precisions.ghs.sd.K10 = unlist(lapply(res.K10.list, FUN = function(s) sd(s$precisions.ghs[,1])))

recalls.sd.K2 = unlist(lapply(res.K2, FUN = function(s) sd(s$recalls[,1])))
recalls.ghs.sd.K2 = unlist(lapply(res.K2, FUN = function(s) sd(s$recalls.ghs[,1])))

recalls.sd.K4 = unlist(lapply(res.K4.list, FUN = function(s) sd(s$recalls[,1])))
recalls.ghs.sd.K4 = unlist(lapply(res.K4.list, FUN = function(s) sd(s$recalls.ghs[,1])))

recalls.sd.K10 = unlist(lapply(res.K10.list, FUN = function(s) sd(s$recalls[,1])))
recalls.ghs.sd.K10 = unlist(lapply(res.K10.list, FUN = function(s) sd(s$recalls.ghs[,1])))


df.K.sd = data.frame(precision=c(precisions.K2, precisions.K2.ghs,precisions.K4,precisions.K4.ghs,precisions.K10,precisions.K10.ghs), 
                  recall=c(recalls.K2, recalls.K2.ghs,recalls.K4,recalls.K4.ghs,recalls.K10,recalls.K10.ghs),
                  sparsity=c(spars.K2, spars.K2.ghs,spars.K4,spars.K4.ghs,spars.K10,spars.K10.ghs),
                  sd.prec = c(precisions.sd.K2, precisions.ghs.sd.K2, precisions.sd.K4, precisions.ghs.sd.K4, precisions.sd.K10, precisions.ghs.sd.K10),
                  sd.rec = c(recalls.sd.K2, recalls.ghs.sd.K2, recalls.sd.K4, recalls.ghs.sd.K4, recalls.sd.K10, recalls.ghs.sd.K10),
                  method = factor(c(rep('jointGHS', length(precisions.K2)),rep('fastGHS',length(precisions.K2.ghs)),
                                    rep('jointGHS', length(precisions.K4)),rep('fastGHS',length(precisions.K4.ghs)),
                                    rep('jointGHS', length(precisions.K10)),rep('fastGHS',length(precisions.K10.ghs))), levels=c('jointGHS','fastGHS')),
                  K = factor(c(rep('K=2', length(precisions.K2)),rep('K=2',length(precisions.K2.ghs)),
                               rep('K=4', length(precisions.K4)),rep('K=4',length(precisions.K4.ghs)),
                               rep('K=10', length(precisions.K10)),rep('K=10',length(precisions.K10.ghs))),levels=c('K=2','K=4','K=10')),
                  disagreement=c(rep(perc.disagreement,6)))

g.prec.sd=ggplot2::ggplot(df.K.sd, aes(y=precision, x=disagreement, group=method, colour=method))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_line(aes(colour=method, linetype=method), size=1)+ scale_color_manual(values=c("darkolivegreen", "darkolivegreen3"))+
  labs(x='Disagreement %')+ scale_x_continuous(breaks = perc.disagreement) + facet_wrap(~K) + #+ scale_y_continuous(breaks = seq(0,60,by=10))
  geom_errorbar(aes(ymin=precision-sd.prec, ymax=precision+sd.prec))+theme_bw()+ theme(legend.position = "none",text = element_text(size = 15))

g.rec.sd = ggplot2::ggplot(df.K.sd, aes(y=recall, x=disagreement, group=method,colour=method))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_line(aes(colour=method, linetype=method), size=1)+ scale_color_manual(values=c("darkolivegreen", "darkolivegreen3"))+
  labs(x='Disagreement %')+ scale_x_continuous(breaks = perc.disagreement) + facet_wrap(~K)+ #+ scale_y_continuous(breaks = seq(0,60,by=10))
  geom_errorbar(aes(ymin=recall-sd.rec, ymax=recall+sd.rec))+ theme_bw()+theme(legend.position = "none",text = element_text(size = 15))

pdf(file='plots/joint_vs_single_SD.pdf',10,8)
gridExtra::grid.arrange(g.prec.sd,g.rec.sd,legend,ncol=1,heights=c(6,6,1))
dev.off()

