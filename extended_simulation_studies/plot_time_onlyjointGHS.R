rm(list=ls())
library(ggplot2)
load("extended_simulation_studies/data/time_simulations_extended_1_joint.Rdata")
time.1 = time.res.large
load("extended_simulation_studies/data/time_simulations_extended_2_joint.Rdata")
time.2 = time.res.large
load("extended_simulation_studies/data/time_simulations_extended_3_joint.Rdata")
time.3 = time.res.large

p=c(600,700,800,900,1000)

df.time.large = data.frame(time=c(time.1,time.2,time.3), p=p)


#pdf(file='extended_simulation_studies/plots/time_jointGHS.pdf',10,8)
pdf(file='extended_simulation_studies/plots/time_jointGHS.pdf',5,5)
ggplot2::ggplot(df.time.large, aes(y=time, x=p))+ labs(title=" ")+theme_bw()+theme(legend.position = 'none', text = element_text(size = 15))+
  geom_line(colour="azure4", linewidth=1)+ 
  labs(y="CPU time (s)")+ scale_x_continuous(breaks = p) +scale_y_continuous(trans='log10')
dev.off()
