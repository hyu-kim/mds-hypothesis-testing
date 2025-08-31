library(superMDS)
library(MASS)
library(parallel)
library(vegan)

#### Import
v_alpha = (0:5)/5

for(rep in 1:3){
  print(sprintf('replicate: %d', rep))
  sim_data <- list()
  sim_data$data <- readRDS(sprintf('result/ScalingStudy/sim_rev_%g/sim_rev_%g-N200-data.Rds', rep, rep))
  sim_data$Y <- read.csv(sprintf('result/ScalingStudy/sim_rev-N200-Y.csv'))
  sim_data$distmat <- as.matrix(vegdist(t(sim_data$data), method="bray"))
  y <- ifelse(sim_data$Y == 1, 1, 2)
  for(alpha in v_alpha){
    res <- TrainSuperMDS(d = sim_data$distmat, y = y, alpha = alpha)
    write.csv(res$z, sprintf('result/HyperparameterStudy/SMDS/sim_rev_%d-smds-%.2f-Z.csv',rep, alpha), row.names=FALSE)
    print('step done')
  }
}