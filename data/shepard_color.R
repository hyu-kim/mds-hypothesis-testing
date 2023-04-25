# source("mm.R")
library(superMDS)
library(scales)
library(grid)
library(gridExtra)
library(parallel)

load('result/sim_etc.RData')

## Simulated superMDS
print(Sys.time())
sim_res_smds <- mclapply(1:4, function(i){
  x <- c(0, 0.3, 0.5, 0.7)[i]
  return(TrainSuperMDS(d = sim_data$distmat, y = sim_data$data$Y, alpha = x/(1+x)))
  print(Sys.time())
}, mc.cores = 4)
names(sim_res_smds) <- c("lambda0", "lambda0.3", "lambda0.5", "lambda0.7")

## Shepard with color
pdf("result/shepard_sim_colored.pdf", width = 7, height = 8)
par(mfrow = c(2, 2))

plot(sim_data$dist, dist(sim_res$mds),
     xlab = "original", ylab = "configuration",
     xaxt="n", yaxt="n",
     main = "proposed, lambda = 0",
     col = alpha("#0944BB", 0.2), 
     pch=16
     )
axis(side=1, at=seq(0,10,5), labels = FALSE)
axis(side=2, at=seq(0,10,5), labels = FALSE)

plot(sim_data$dist, dist(sim_res$proposed$lambda0.7$z),
     xlab = "original", ylab = "configuration",
     xaxt="n", yaxt="n",
     main = "proposed, lambda = 0.7", 
     col = alpha("#0944BB", 0.2), 
     pch=16
     )
axis(side=1, at=seq(0,10,5), labels = FALSE)
axis(side=2, at=seq(0,10,5), labels = FALSE)

plot(sim_data$dist, dist(sim_res_smds$lambda0$z),
     xlab = "original", ylab = "configuration",
     xaxt="n", yaxt="n",
     main = "SMDS, lambda = 0", 
     col = alpha("#DB4915", 0.2), 
     pch=16
     )
axis(side=1, at=seq(0,10,5), labels = FALSE)
axis(side=2, at=seq(0,10,5), labels = FALSE)

plot(sim_data$dist, dist(sim_res_smds$lambda0.7$z),
     xlab = "original", ylab = "configuration",
     xaxt="n", yaxt="n",
     main = "SMDS, lambda = 0.7", 
     col = alpha("#DB4915", 0.15), 
     pch=16
     )
axis(side=1, at=seq(0,10,5), labels = FALSE)
axis(side=2, at=seq(0,10,5), labels = FALSE)

dev.off()