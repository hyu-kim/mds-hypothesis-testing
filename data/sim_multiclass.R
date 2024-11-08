#3. Simulation
library(MASS)
library(parallel)
library(pbmcapply)
source('mm.R')

set.seed(150)
sim_data <- list(
  data = rbind(
    mvrnorm(50, c(0,0,0,0), 
            Sigma = matrix(c(5,0,0,0, 0,5,0,0, 0,0,1,0, 0,0,0,1), nrow=4)),
    mvrnorm(50, c(0,0,2,0), 
            Sigma = matrix(c(5,0,0,0, 0,5,0,0, 0,0,1,0, 0,0,0,1), nrow=4)),
    mvrnorm(50, c(0,0,1,sqrt(3)), 
            Sigma = matrix(c(5,0,0,0, 0,5,0,0, 0,0,1,0, 0,0,0,1), nrow=4))
    )
)
sim_data$data <- data.frame(sim_data$data)
sim_data$Y <- c(rep(1, 50), rep(2, 50), rep(3, 50))

sim_data$distmat <- as.matrix(dist(sim_data$data[, 1:4]))

z0 <- cmdscale(dist(sim_data$data[,1:4]), k = 2)

write.csv(sim_data$data, 'result/Multiclass/sim4d_2-data.csv', row.names=FALSE)
write.csv(sim_data$Y, 'result/Multiclass/sim4d_2-Y.csv', row.names=FALSE)
write.csv(z0, 'result/Multiclass/sim4d_2-mds-Z.csv', row.names=FALSE)

res <- mm_cmds(nit=100, lambda=0.5, z0=z0, D=sim_data$distmat, y=sim_data$Y, dataset = 'sim4d_2')

