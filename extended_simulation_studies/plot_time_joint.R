rm(list=ls())
library(ggplot2)
source('simulation_functions/help_functions.R')
source('simulation_functions/perform_time_simulation.R')

# Plot results from time simulations

load("extended_simulation_studies/data/time_simulations_joint_K2.Rdata")
load("extended_simulation_studies/data/time_simulations_joint_K3.Rdata")
load("extended_simulation_studies/data/time_simulations_joint_K4.Rdata")
p=c(20,50,100,150,200,250)
n.p = length(p)



# All joint methods

df.time.K2 = data.frame(time=c(time.res.joint.K2$times.jointGHS, time.res.joint.K2$times.JGL, time.res.joint.K2$times.SSJGL, time.res.joint.K2$times.GemBag), p=rep(p,2),
                        method=c(rep('jointGHS',n.p),rep('JGL',n.p),rep('SSJGL',n.p),rep('GemBag',n.p) ))

df.time.K3 = data.frame(time=c(time.res.joint.K3$times.jointGHS, time.res.joint.K3$times.JGL, time.res.joint.K3$times.SSJGL, time.res.joint.K3$times.GemBag), p=rep(p,2),
                        method=c(rep('jointGHS',n.p),rep('JGL',n.p),rep('SSJGL',n.p),rep('GemBag',n.p) ))

df.time.K4 = data.frame(time=c(time.res.joint.K4$times.jointGHS, time.res.joint.K4$times.JGL, time.res.joint.K4$times.SSJGL, time.res.joint.K4$times.GemBag), p=rep(p,2),
                        method=c(rep('jointGHS',n.p),rep('JGL',n.p),rep('SSJGL',n.p),rep('GemBag',n.p) ))


df.time = rbind(df.time.K2,df.time.K3,df.time.K4)
df.time$K=c(rep('K=2', n.p*4),rep('K=3', n.p*4),rep('K=4', n.p*4))

pdf(file='extended_simulation_studies/plots/time_joint.pdf',8,4)
ggplot2::ggplot(df.time, aes(y=time, x=p, group=method))+ labs(title=" ")+theme_bw()+theme(plot.title = element_text(hjust = 0.5),text = element_text(size = 15))+
  geom_line(aes(colour=method, linetype=method), size=1)+ scale_color_manual(values=c("azure4", "navyblue","firebrick3","palegreen4"))+
  labs(y="CPU time (s)")+ scale_x_continuous(breaks = p) +facet_wrap(~K) #+ scale_y_continuous(breaks = seq(0,60,by=10))
dev.off()

pdf(file='extended_simulation_studies/plots/time_joint_log.pdf',8,4)
ggplot2::ggplot(df.time, aes(y=time, x=p, group=method))+ labs(title=" ")+theme_bw()+theme(plot.title = element_text(hjust = 0.5),text = element_text(size = 15))+
  geom_line(aes(colour=method, linetype=method), linewidth=1)+ scale_color_manual(values=c("azure4", "navyblue","firebrick3","palegreen4"))+
  labs(y="CPU time (s)")+ scale_x_continuous(breaks = p)+scale_y_continuous(trans='log10',breaks = c(0,5, 50,500,5000,50000))+facet_wrap(~K)
dev.off()









