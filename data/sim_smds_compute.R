library(superMDS)
library(MASS)
library(parallel)

#### Import
v_alpha = (0:5)/5

for(rep in 1:3){
  print(sprintf('replicate: %d', rep))
  sim_data <- list()
  sim_data$data <- read.csv(paste('data/sim/sim-data-',rep,'.csv', sep=''))
  sim_data$Y <- read.csv(paste('data/sim/sim-Y-',rep,'.csv', sep=''))
  y <- ifelse(sim_data$Y == 1, 1, 2)
  sim_data$distmat <- as.matrix(dist(sim_data$data[, 1:3]))
  for(alpha in v_alpha){
    res <- TrainSuperMDS(d = sim_data$distmat, y = y, alpha = alpha)
    write.csv(res$z, sprintf('result/smds_sim-%.2f-Z-%d.csv', alpha, rep), row.names=FALSE)
    print('step done')
  }
}