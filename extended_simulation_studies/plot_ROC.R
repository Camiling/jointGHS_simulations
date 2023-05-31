rm(list=ls())
library(ggplot2)
library(dplyr)
source('simulation_functions/help_functions.R')
load("extended_simulation_studies/data/jointGHS_simulations_jointGHS_ROC.Rdata")
load("extended_simulation_studies/data/jointGHS_simulations_JGL_ROC.Rdata")
load("extended_simulation_studies/data/jointGHS_simulations_SSJGL_ROC.Rdata")
load("extended_simulation_studies/data/jointGHS_simulations_GemBag_ROC.Rdata")

# We limit ourselves to disagreement levels 0%, 60% and 100% (cases 1, 4 and 6). 

# PRECISION-RECALL CURVES ---------------------

n.jointGHS = length(res.ROC.jointGHS[[1]]$mean.precisions[2,])
n.SSJGL = length(res.ROC.SSJGL[[1]]$mean.precisions.ssjgl[2,])
n.JGL = length(res.ROC.JGL[[1]]$mean.precisions.jgl[2,])
n.gembag = length(res.ROC.gembag[[1]]$mean.precisions.gembag[2,])

# 0% disagreement
plots.0 = list()
for(j in 1:2){
  df.0 = data.frame(precision = c(res.ROC.jointGHS[[1]]$mean.precisions[j,],res.ROC.SSJGL[[1]]$mean.precisions.ssjgl[j,],res.ROC.JGL[[1]]$mean.precisions.jgl[j,],res.ROC.gembag[[1]]$mean.precisions.gembag[j,]),
                    recall = c(res.ROC.jointGHS[[1]]$mean.recalls[j,],res.ROC.SSJGL[[1]]$mean.recalls.ssjgl[j,],res.ROC.JGL[[1]]$mean.recalls.jgl[j,],res.ROC.gembag[[1]]$mean.recalls.gembag[j,]),
                    method = c(rep('jointGHS', n.jointGHS), rep('SSJGL', n.SSJGL),rep('JGL', n.JGL),rep('GemBag', n.gembag)))
  df.0 = df.0 %>% filter(recall <= max(res.ROC.jointGHS[[1]]$mean.recalls[j,])+0.01)
  plots.0[[j]] = ggplot(df.0, aes(x=recall, y=precision, color=method, linetype=method))+geom_line(data= df.0 %>% filter(method %in% c('jointGHS','SSJGL','JGL')),linewidth=1)+
    geom_point(data= df.0 %>% filter(method %in% c('GemBag')), size=2)+ theme_bw()+ scale_color_manual(values=c('thistle3','lightskyblue3','dodgerblue','darkkhaki'))+
    ggtitle('0% disagreement')+theme(legend.position = 'none',plot.title = element_text(hjust = 0.5),text = element_text(size = 15))+ ylim(0,1)+ xlim(0,0.3)
}

# 40% disagreement
plots.40 = list()
for(j in 1:2){
  df.40 = data.frame(precision = c(res.ROC.jointGHS[[3]]$mean.precisions[j,],res.ROC.SSJGL[[3]]$mean.precisions.ssjgl[j,],res.ROC.JGL[[3]]$mean.precisions.jgl[j,],res.ROC.gembag[[3]]$mean.precisions.gembag[j,]),
                     recall = c(res.ROC.jointGHS[[3]]$mean.recalls[j,],res.ROC.SSJGL[[3]]$mean.recalls.ssjgl[j,],res.ROC.JGL[[3]]$mean.recalls.jgl[j,],res.ROC.gembag[[3]]$mean.recalls.gembag[j,]),
                     method = c(rep('jointGHS', n.jointGHS), rep('SSJGL', n.SSJGL),rep('JGL', n.JGL),rep('GemBag', n.gembag)))
  df.40 = df.40 %>% filter(recall <= max(res.ROC.jointGHS[[3]]$mean.recalls[j,])+0.01)
  plots.40[[j]] = ggplot(df.40, aes(x=recall, y=precision, color=method, linetype=method))+geom_line(data= df.40 %>% filter(method %in% c('jointGHS','SSJGL','JGL')),linewidth=1)+
    geom_point(data= df.40 %>% filter(method %in% c('GemBag')), size=2)+ theme_bw()+ scale_color_manual(values=c('thistle3','lightskyblue3','dodgerblue','darkkhaki'))+
    ggtitle('40% disagreement')+theme(legend.position = 'none',plot.title = element_text(hjust = 0.5),text = element_text(size = 15))+ ylim(0,1)+ xlim(0,0.3)
}

# 60% disagreement
plots.60 = list()
for(j in 1:2){
  df.60 = data.frame(precision = c(res.ROC.jointGHS[[4]]$mean.precisions[j,],res.ROC.SSJGL[[4]]$mean.precisions.ssjgl[j,],res.ROC.JGL[[4]]$mean.precisions.jgl[j,],res.ROC.gembag[[4]]$mean.precisions.gembag[j,]),
                     recall = c(res.ROC.jointGHS[[4]]$mean.recalls[j,],res.ROC.SSJGL[[4]]$mean.recalls.ssjgl[j,],res.ROC.JGL[[4]]$mean.recalls.jgl[j,],res.ROC.gembag[[4]]$mean.recalls.gembag[j,]),
                     method = c(rep('jointGHS', n.jointGHS), rep('SSJGL', n.SSJGL),rep('JGL', n.JGL),rep('GemBag', n.gembag)))
  df.60 = df.60 %>% filter(recall <= max(res.ROC.jointGHS[[4]]$mean.recalls[j,])+0.01)
  plots.60[[j]] = ggplot(df.60, aes(x=recall, y=precision, color=method, linetype=method))+geom_line(data= df.60 %>% filter(method %in% c('jointGHS','SSJGL','JGL')),linewidth=1)+
    geom_point(data= df.60 %>% filter(method %in% c('GemBag')), size=2)+ theme_bw()+ scale_color_manual(values=c('thistle3','lightskyblue3','dodgerblue','darkkhaki'))+
    ggtitle('60% disagreement')+theme(legend.position = 'none',plot.title = element_text(hjust = 0.5),text = element_text(size = 15))+ ylim(0,1)+ xlim(0,0.3)
}

# 100% disagreement
plots.100 = list()
for(j in 1:2){
  df.100 = data.frame(precision = c(res.ROC.jointGHS[[6]]$mean.precisions[j,],res.ROC.SSJGL[[6]]$mean.precisions.ssjgl[j,],res.ROC.JGL[[6]]$mean.precisions.jgl[j,],res.ROC.gembag[[6]]$mean.precisions.gembag[j,]),
                      recall = c(res.ROC.jointGHS[[6]]$mean.recalls[j,],res.ROC.SSJGL[[6]]$mean.recalls.ssjgl[j,],res.ROC.JGL[[6]]$mean.recalls.jgl[j,],res.ROC.gembag[[6]]$mean.recalls.gembag[j,]),
                      method = c(rep('jointGHS', n.jointGHS), rep('SSJGL', n.SSJGL),rep('JGL', n.JGL),rep('GemBag', n.gembag)))
  df.100 = df.100 %>% filter(recall <= max(res.ROC.jointGHS[[6]]$mean.recalls[j,])+0.01)
  plots.100[[j]] = ggplot(df.100, aes(x=recall, y=precision, color=method, linetype=method))+geom_line(data= df.100 %>% filter(method %in% c('jointGHS','SSJGL','JGL')),linewidth=1)+
    geom_point(data= df.100 %>% filter(method %in% c('GemBag')), size=2)+ theme_bw()+ scale_color_manual(values=c('thistle3','lightskyblue3','dodgerblue','darkkhaki'))+
    ggtitle('100% disagreement')+theme(plot.title = element_text(hjust = 0.5),text = element_text(size = 15), legend.position = 'bottom')+ ylim(0,1)+ xlim(0,0.3)
}

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend=get_legend(plots.100[[2]])

layout.mat = matrix(c(1,2,
                      3,4,
                      5,5),ncol=2,byrow=T)

pdf(file='extended_simulation_studies/plots/PRC_ROC_curves_g2.pdf',8,6)
#gridExtra::grid.arrange(plots.0[[2]], plots.40[[2]],legend,plots.60[[2]],plots.100[[2]]+theme(legend.position = 'none',plot.title = element_text(hjust = 0.5),text = element_text(size = 15)),ncol=3,widths=c(1,1,0.3))
gridExtra::grid.arrange(plots.0[[2]], plots.40[[2]],plots.60[[2]],plots.100[[2]]+theme(legend.position = 'none',plot.title = element_text(hjust = 0.5),text = element_text(size = 15)),legend,
                        layout_matrix=layout.mat,heights=c(1,1,0.2))
dev.off()



# COMPUTE AUCPRC ------------------------------

get.AUPRC.f = function(rec,prec,rec.cut=0.01){
  prec=prec[order(rec)]
  rec=rec[order(rec)]
  rec = rec[which(rec<=rec.cut)]
  n.points = length(rec)
  auc = sum((rec[2:n.points]-rec[1:(n.points-1)])*prec[2:n.points])
  # Extend end values so we include 0 and f
  auc = auc + (rec.cut-rec[n.points])*prec[n.points] + (rec[1]-0)*prec[1]
  return(auc)
}

# We rescale the scores according to f

# 0% disagreement
#f=max(res.ROC.jointGHS[[1]]$mean.recalls[2,])
f=0.3
f
round(get.AUPRC.f(res.ROC.jointGHS[[1]]$mean.recalls[2,],res.ROC.jointGHS[[1]]$mean.precisions[2,],f), 3)
# 0.3
round(get.AUPRC.f(res.ROC.JGL[[1]]$mean.recalls.jgl[2,],res.ROC.JGL[[1]]$mean.precisions.jgl[2,],f), 3)
# 0.295
round(get.AUPRC.f(res.ROC.SSJGL[[1]]$mean.recalls.ssjgl[2,],res.ROC.SSJGL[[1]]$mean.precisions.ssjgl[2,],f), 3)
# 0.292

# 40% disagreement
#f=max(res.ROC.jointGHS[[3]]$mean.recalls[2,])
f=0.3
f
round(get.AUPRC.f(res.ROC.jointGHS[[3]]$mean.recalls[2,],res.ROC.jointGHS[[3]]$mean.precisions[2,],f), 3)
# 0.299
round(get.AUPRC.f(res.ROC.JGL[[3]]$mean.recalls.jgl[2,],res.ROC.JGL[[3]]$mean.precisions.jgl[2,],f), 3)
# 0.287
round(get.AUPRC.f(res.ROC.SSJGL[[3]]$mean.recalls.ssjgl[2,],res.ROC.SSJGL[[3]]$mean.precisions.ssjgl[2,],f), 3)
# 0.286

# 60% disagreement
#f=max(res.ROC.jointGHS[[4]]$mean.recalls[2,])
f=0.3
f
round(get.AUPRC.f(res.ROC.jointGHS[[4]]$mean.recalls[2,],res.ROC.jointGHS[[4]]$mean.precisions[2,],f), 3)
# 0.299
round(get.AUPRC.f(res.ROC.JGL[[4]]$mean.recalls.jgl[2,],res.ROC.JGL[[4]]$mean.precisions.jgl[2,],f), 3)
# 0.285
round(get.AUPRC.f(res.ROC.SSJGL[[4]]$mean.recalls.ssjgl[2,],res.ROC.SSJGL[[4]]$mean.precisions.ssjgl[2,],f), 3)
# 0.285

# 100% disagreement
#f=max(res.ROC.jointGHS[[6]]$mean.recalls[2,])
f=0.3
f
round(get.AUPRC.f(res.ROC.jointGHS[[6]]$mean.recalls[2,],res.ROC.jointGHS[[6]]$mean.precisions[2,],f), 3)
# 0.299
round(get.AUPRC.f(res.ROC.JGL[[6]]$mean.recalls.jgl[2,],res.ROC.JGL[[6]]$mean.precisions.jgl[2,],f), 3)
# 0.263
round(get.AUPRC.f(res.ROC.SSJGL[[6]]$mean.recalls.ssjgl[2,],res.ROC.SSJGL[[6]]$mean.precisions.ssjgl[2,],f), 3)
# 0.267







