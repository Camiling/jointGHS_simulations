source('simulation_functions/perform_jointGHS_simulation.R')
source('simulation_functions/help_functions.R')
load("data/jointGHS_simulations_K2.Rdata")

# 100 simulations per case
N = 100
p = 50
K=2
n.vals = c(50,80)
fracs.disagreement = c(0,0.2,0.4,0.6,0.8,1)
n.points = length(fracs.disagreement)
perc.disagreement = 100*fracs.disagreement


# Plot results

# First network
precisions.K2 = unlist(lapply(res.K2, FUN = function(s) mean(s$mean.precisions[1])))
precisions.K2.ssjgl = unlist(lapply(res.K2, FUN = function(s) mean(s$mean.precisions.ssjgl[1])))
precisions.K2.jgl = unlist(lapply(res.K2, FUN = function(s) mean(s$mean.precisions.jgl[1])))

recalls.K2 = unlist(lapply(res.K2, FUN = function(s) mean(s$mean.recalls[1])))
recalls.K2.ssjgl = unlist(lapply(res.K2, FUN = function(s) mean(s$mean.recalls.ssjgl[1])))
recalls.K2.jgl = unlist(lapply(res.K2, FUN = function(s) mean(s$mean.recalls.jgl[1])))

spars.K2 = unlist(lapply(res.K2, FUN = function(s) mean(s$mean.opt.sparsities[1])))
spars.K2.ssjgl = unlist(lapply(res.K2, FUN = function(s) mean(s$mean.opt.sparsities.ssjgl[1])))
spars.K2.jgl = unlist(lapply(res.K2, FUN = function(s) mean(s$mean.opt.sparsities.jgl[1])))

# Standard deviations
precisions.sd.K2 = unlist(lapply(res.K2, FUN = function(s) sd(s$precisions[,1])))
precisions.ssjgl.sd.K2 = unlist(lapply(res.K2, FUN = function(s) sd(s$precisions.ssjgl[,1])))
precisions.jgl.sd.K2 = unlist(lapply(res.K2, FUN = function(s) sd(s$precisions.jgl[,1])))

recalls.sd.K2 = unlist(lapply(res.K2, FUN = function(s) sd(s$recalls[,1])))
recalls.ssjgl.sd.K2 = unlist(lapply(res.K2, FUN = function(s) sd(s$recalls.ssjgl[,1])))
recalls.jgl.sd.K2 = unlist(lapply(res.K2, FUN = function(s) sd(s$recalls.jgl[,1])))

# Second network
precisions.K2.2 = unlist(lapply(res.K2, FUN = function(s) mean(s$mean.precisions[2])))
precisions.K2.2.ssjgl = unlist(lapply(res.K2, FUN = function(s) mean(s$mean.precisions.ssjgl[2])))
precisions.K2.2.jgl = unlist(lapply(res.K2, FUN = function(s) mean(s$mean.precisions.jgl[2])))

recalls.K2.2 = unlist(lapply(res.K2, FUN = function(s) mean(s$mean.recalls[2])))
recalls.K2.2.ssjgl = unlist(lapply(res.K2, FUN = function(s) mean(s$mean.recalls.ssjgl[2])))
recalls.K2.2.jgl = unlist(lapply(res.K2, FUN = function(s) mean(s$mean.recalls.jgl[2])))

spars.K2.2 = unlist(lapply(res.K2, FUN = function(s) mean(s$mean.opt.sparsities[2])))
spars.K2.2.ssjgl = unlist(lapply(res.K2, FUN = function(s) mean(s$mean.opt.sparsities.ssjgl[2])))
spars.K2.2.jgl = unlist(lapply(res.K2, FUN = function(s) mean(s$mean.opt.sparsities.jgl[2])))

# Standard deviations
precisions.sd.K2.2 = unlist(lapply(res.K2, FUN = function(s) sd(s$precisions[,2])))
precisions.ssjgl.sd.K2.2 = unlist(lapply(res.K2, FUN = function(s) sd(s$precisions.ssjgl[,2])))
precisions.jgl.sd.K2.2 = unlist(lapply(res.K2, FUN = function(s) sd(s$precisions.jgl[,2])))

recalls.sd.K2.2 = unlist(lapply(res.K2, FUN = function(s) sd(s$recalls[,2])))
recalls.ssjgl.sd.K2.2 = unlist(lapply(res.K2, FUN = function(s) sd(s$recalls.ssjgl[,2])))
recalls.jgl.sd.K2.2 = unlist(lapply(res.K2, FUN = function(s) sd(s$recalls.jgl[,2])))


df.joint = data.frame(precision=c(precisions.K2, precisions.K2.ssjgl, precisions.K2.jgl, precisions.K2.2, precisions.K2.2.ssjgl, precisions.K2.2.jgl),
                      recall = c(recalls.K2, recalls.K2.ssjgl, recalls.K2.jgl, recalls.K2.2, recalls.K2.2.ssjgl, recalls.K2.2.jgl),
                      sparsity = c(spars.K2, spars.K2.ssjgl, spars.K2.jgl, spars.K2.2, spars.K2.2.ssjgl, spars.K2.2.jgl),
                      sd.prec = c(precisions.sd.K2, precisions.ssjgl.sd.K2, precisions.jgl.sd.K2, precisions.sd.K2.2, precisions.ssjgl.sd.K2.2, precisions.jgl.sd.K2.2),
                      sd.rec = c(recalls.sd.K2, recalls.ssjgl.sd.K2, recalls.jgl.sd.K2, recalls.sd.K2.2, recalls.ssjgl.sd.K2.2, recalls.jgl.sd.K2.2),
                      method = factor(rep(c(rep('jointGHS', n.points), rep('SSJGL', n.points), rep('JGL', n.points)), 2),levels = c('jointGHS', 'SSJGL', 'JGL')), 
                      disagreement=c(rep(perc.disagreement,6)), 
                      graph = factor(c(rep('Graph 1', n.points*3),rep('Graph 2', n.points*3)), levels=c('Graph 1', 'Graph 2')))


g.prec = ggplot2::ggplot(df.joint, aes(y=precision, x=disagreement, group=method, colour=method))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_line(aes(colour=method, linetype=method), size=1)+ scale_color_manual(values=c("darkolivegreen", "darkblue", "lightblue"))+
  labs(x='Disagreement %')+ scale_x_continuous(breaks = perc.disagreement) + facet_wrap(~graph) + #+ scale_y_continuous(breaks = seq(0,60,by=10))
   theme(legend.position = "bottom")

g.rec = ggplot2::ggplot(df.joint, aes(y=recall, x=disagreement, group=method,colour=method))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_line(aes(colour=method, linetype=method), size=1)+ scale_color_manual(values=c("darkolivegreen", "darkblue", "lightblue"))+
  labs(x='Disagreement %')+ scale_x_continuous(breaks = perc.disagreement) + facet_wrap(~graph)


g.prec.sd = ggplot2::ggplot(df.joint, aes(y=precision, x=disagreement, group=method, colour=method))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_line(aes(colour=method, linetype=method), size=1)+ scale_color_manual(values=c("darkolivegreen", "darkblue", "lightblue"))+
  labs(x='Disagreement %')+ scale_x_continuous(breaks = perc.disagreement) + facet_wrap(~graph) + #+ scale_y_continuous(breaks = seq(0,60,by=10))
  geom_errorbar(aes(ymin=precision-sd.prec, ymax=precision+sd.prec)) + theme(legend.position = "bottom")

g.rec.sd = ggplot2::ggplot(df.joint, aes(y=recall, x=disagreement, group=method,colour=method))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_line(aes(colour=method, linetype=method), size=1)+ scale_color_manual(values=c("darkolivegreen", "darkblue", "lightblue"))+
  labs(x='Disagreement %')+ scale_x_continuous(breaks = perc.disagreement) + facet_wrap(~graph)+ #+ scale_y_continuous(breaks = seq(0,60,by=10))
  geom_errorbar(aes(ymin=recall-sd.rec, ymax=recall+sd.rec))


tmp = ggplot2::ggplot_gtable(ggplot2::ggplot_build(g.prec+theme(text = element_text(size = 15))))
leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
legend = tmp$grobs[[leg]]
g.rec.sd = g.rec.sd + theme_bw()+theme(legend.position = "none",text = element_text(size = 15))
g.prec.sd = g.prec.sd + theme_bw()+theme(legend.position = "none",text = element_text(size = 15))
g.rec = g.rec + theme_bw()+theme(legend.position = "none",text = element_text(size = 15))
g.prec = g.prec + theme_bw()+theme(legend.position = "none",text = element_text(size = 15))


pdf(file='plots/joint_comparison_all_SD.pdf',10,8)
gridExtra::grid.arrange(g.prec.sd,g.rec.sd,legend,ncol=1,heights=c(6,6,1))
dev.off()

pdf(file='plots/joint_comparison.pdf',10,8)
gridExtra::grid.arrange(g.prec,g.rec,legend,ncol=1,heights=c(6,6,1))
dev.off()

# Just the precision

pdf(file='plots/joint_comparison_all_prec.pdf',10,7)
gridExtra::grid.arrange(g.prec.sd,legend,ncol=1,heights=c(7,1))
dev.off()

# plot without JGL

g.prec.nojgl = ggplot2::ggplot(df.joint[df.joint$method!='JGL',], aes(y=precision, x=disagreement, group=method, colour=method))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_line(aes(colour=method, linetype=method), size=1)+ scale_color_manual(values=c("darkolivegreen", "darkblue"))+
  labs(x='Disagreement %')+ scale_x_continuous(breaks = perc.disagreement) + facet_wrap(~graph) + #+ scale_y_continuous(breaks = seq(0,60,by=10))
  theme(legend.position = "bottom")

g.rec.nojgl  = ggplot2::ggplot(df.joint[df.joint$method!='JGL',], aes(y=recall, x=disagreement, group=method,colour=method))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_line(aes(colour=method, linetype=method), size=1)+ scale_color_manual(values=c("darkolivegreen", "darkblue"))+
  labs(x='Disagreement %')+ scale_x_continuous(breaks = perc.disagreement) + facet_wrap(~graph)


g.prec.sd.nojgl  = ggplot2::ggplot(df.joint[df.joint$method!='JGL',], aes(y=precision, x=disagreement, group=method, colour=method))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_line(aes(colour=method, linetype=method), size=1)+ scale_color_manual(values=c("darkolivegreen", "darkblue"))+
  labs(x='Disagreement %')+ scale_x_continuous(breaks = perc.disagreement) + facet_wrap(~graph) + #+ scale_y_continuous(breaks = seq(0,60,by=10))
  geom_errorbar(aes(ymin=precision-sd.prec, ymax=precision+sd.prec)) + theme(legend.position = "bottom")

g.rec.sd.nojgl  = ggplot2::ggplot(df.joint[df.joint$method!='JGL',], aes(y=recall, x=disagreement, group=method,colour=method))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_line(aes(colour=method, linetype=method), size=1)+ scale_color_manual(values=c("darkolivegreen", "darkblue"))+
  labs(x='Disagreement %')+ scale_x_continuous(breaks = perc.disagreement) + facet_wrap(~graph)+ #+ scale_y_continuous(breaks = seq(0,60,by=10))
  geom_errorbar(aes(ymin=recall-sd.rec, ymax=recall+sd.rec))

tmp = ggplot2::ggplot_gtable(ggplot2::ggplot_build(g.prec.nojgl+theme(text = element_text(size = 15))))
leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
legend = tmp$grobs[[leg]]
g.rec.sd.nojgl = g.rec.sd.nojgl + theme_bw()+theme(legend.position = "none",text = element_text(size = 15))
g.prec.sd.nojgl = g.prec.sd.nojgl + theme_bw()+theme(legend.position = "none",text = element_text(size = 15))
g.rec.nojgl = g.rec.nojgl + theme_bw()+theme(legend.position = "none",text = element_text(size = 15))
g.prec.nojgl = g.prec.nojgl + theme_bw()+theme(legend.position = "none",text = element_text(size = 15))


pdf(file='plots/joint_comparison_small_SD.pdf',10,8)
gridExtra::grid.arrange(g.prec.sd.nojgl,g.rec.sd.nojgl,legend,ncol=1,heights=c(6,6,1))
dev.off()

pdf(file='plots/joint_comparison_small.pdf',10,8)
gridExtra::grid.arrange(g.prec.nojgl,g.rec.nojgl,legend,ncol=1,heights=c(6,6,1))
dev.off()








