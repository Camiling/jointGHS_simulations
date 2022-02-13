source('simulation_functions/perform_fastGHS_simulation.R')
source('simulation_functions/help_functions.R')
source('simulation_functions/perform_time_simulation.R')

# Plot results from time simulations

load("data/time_simulations_small.Rdata")
load("data/time_simulations_large.Rdata")
p=c(10,20,30,40,50,60,70,80,90)

# GHS vs fastGHS

df.time = data.frame(time=c(time.res[[1]], time.res[[2]]), p=rep(p,2),
                     method=c(rep('fastGHS',length(time.res[[1]])),rep('GHS',length(time.res[[2]])) ))


pdf(file='plots/time_small.pdf',8,8)
ggplot2::ggplot(df.time, aes(y=time, x=p, group=method))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_line(aes(colour=method, linetype=method))+ scale_color_manual(values=c("azure4", "darkgray"))+
  labs(y="CPU time (s)")+ scale_x_continuous(breaks = p)+ scale_y_continuous(breaks = seq(0,60,by=10))
dev.off()

# Only fastGHS, larger p

p=c(20,50,100,200,300,400,500,600)
df.time.large = data.frame(time=c(time.res.large), p=p,
                     method=c(rep('fastGHS',length(time.res.large))))


pdf(file='plots/time_large.pdf',8,8)
ggplot2::ggplot(df.time.large, aes(y=time, x=p, group=method))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_line(aes(colour=method, linetype=method))+ scale_color_manual(values=c("azure4"))+
  labs(y="CPU time (s)")+ scale_x_continuous(breaks = p[3:8]) #+ scale_y_continuous(breaks = seq(0,60,by=10))
dev.off()

pdf(file='plots/time_large_nolabel.pdf',8,8)
ggplot2::ggplot(df.time.large, aes(y=time, x=p, group=method))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_line(aes(colour=method, linetype=method))+ scale_color_manual(values=c("azure4"))+
  labs(y="CPU time (s)")+ scale_x_continuous(breaks = p[3:8]) + theme(legend.position = "none") #+ scale_y_continuous(breaks = seq(0,60,by=10))
dev.off()

# Plot time per edge

df.time.large.per = df.time.large
df.time.large.per$time = df.time.large.per$time/((df.time.large.per$p^2-df.time.large.per$p)/2)

pdf(file='plots/time_large_pernode.pdf',8,8)
ggplot2::ggplot(df.time.large.per, aes(y=time, x=p, group=method))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_line(aes(colour=method, linetype=method))+ scale_color_manual(values=c("azure4"))+
  labs(y="CPU time (s)")+ scale_x_continuous(breaks = p[3:8]) #+ scale_y_continuous(breaks = seq(0,60,by=10))
dev.off()

pdf(file='plots/time_large_pernode_nolabel.pdf',8,8)
ggplot2::ggplot(df.time.large.per, aes(y=time, x=p, group=method))+ labs(title=" ")+theme(plot.title = element_text(hjust = 0.5))+
  geom_line(aes(colour=method, linetype=method))+ scale_color_manual(values=c("azure4"))+
  labs(y="CPU time (s)")+ scale_x_continuous(breaks = p[3:8]) + theme(legend.position = "none")#+ scale_y_continuous(breaks = seq(0,60,by=10))
dev.off()

