#3. Simulation
library(MASS)
library(parallel)
library(pbmcapply)
set.seed(100)
sim_data <- list(
  data = rbind(
    mvrnorm(50, c(0,0,0), Sigma = matrix(c(3,0,0,0,3,0,0,0,1), nrow=3)),
    mvrnorm(50, c(0,0,1), Sigma = matrix(c(3,0,0,0,3,0,0,0,1), nrow=3)))
)
sim_data$data <- data.frame(sim_data$data)
sim_data$Y <- c(rep(1, 50), rep(2, 50))

sim_data$distmat <- as.matrix(dist(sim_data$data[, 1:3]))

z0 <- cmdscale(dist(sim_data$data[,1:3]), k = 2)

write.csv(sim_data$data, 'result/sim-data.csv', row.names=FALSE)
write.csv(sim_data$Y, 'result/sim-Y.csv', row.names=FALSE)

print('exporting done. Begins iteration..')
res <- pbmclapply(1:21, function(i){
  x <- c(0:20)[i]/20
  return(mm_cmds(nit=50, 
                 lambda=x, 
                 z0=z0, 
                 D=sim_data$distmat, 
                 y = sim_data$Y,
                 dataset = 'sim')
         )
}, mc.cores = 128)

# lambda-free
source('mm.R')
print('lambda-free method')
model_loess <- readRDS('result/HyperparameterStudy/Nonparametric/model_loess.Rds')
res <- mm_cmds2(nit=200, z0=z0, D=sim_data$distmat, y = sim_data$Y, 
                dataset = 'sim', model_lambda=model_loess)

write.csv(data.frame(n_iter=200, dataset='sim'), 'result/sim-config.csv', row.names=FALSE)
