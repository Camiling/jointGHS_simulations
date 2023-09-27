source('simulation_functions/perform_fastGHS_simulation.R')
source('simulation_functions/help_functions.R')

# Write print function

nCores = 20 # If using HPC

# Perform simulation study

perform_GHS_sims = TRUE

# 20 simulations per case
N = 100

if(perform_GHS_sims){
  # Case 1: p=50,n=100
  set.seed(12345)
  p = 50
  n = 100
  res.1 = perform_fastGHS_simulation(n, p, N, nCores = nCores, include.glasso = F,include.GHS = F, save_Q=T)
  
  # Case 2: p=50, n=200
  set.seed(1234)
  p = 50
  n = 200
  res.2 = perform_fastGHS_simulation(n, p, N, nCores = nCores, include.glasso = F,include.GHS = F, save_Q=T)
  
  # Case 3: p=100, n=100
  set.seed(1234)
  p = 100
  n = 100
  res.3 = perform_fastGHS_simulation(n, p, N, nCores = nCores, include.glasso = F,include.GHS = F, save_Q=T)
  
  # Case 4: p=100, n=200
  set.seed(1234)
  p = 100
  n = 200
  res.4 = perform_fastGHS_simulation(n, p, N, nCores = nCores, include.glasso = F,include.GHS = F, save_Q=T)
  
  # Save results
  res.fast = list(res.1, res.2, res.3, res.4)
  save(res.fast, file="extended_simulation_studies/data/fastGHS_simulations_randominit.Rdata")
  load("extended_simulation_studies/data/fastGHS_simulations_randominit.Rdata")
  
}


# Plot results

df.init = data.frame(precision=unlist(lapply(res.fast,FUN=function(s) s$precisions.fastghs)),
                     recall = unlist(lapply(res.fast,FUN=function(s) s$recalls.fastghs)),
                     setting = c(rep('p=50, n=100',N), rep('p=50, n=200',N), rep('p=100, n=100',N),rep('p=100, n=200',N))
                     )
  
pdf(file='extended_simulation_studies/plots/plot_init.pdf',6,6)
ggplot2::ggplot(df.init, aes(x=recall, y=precision, group=setting))+geom_point(size=2,colour='darkgray')+ facet_wrap(~setting)+xlim(0,1)+ylim(0,1)+
                theme_bw() + theme(text = element_text(size = 15), legend.position = 'none')
dev.off()
  
  
  
  
  

  
  


