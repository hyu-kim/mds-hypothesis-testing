library(superMDS)

#### Simulation
set.seed(10101)
library(MASS)
x <- mvrnorm(100, rep(0, 10), diag(1, 10, 10))
y_lab <- c(rep(1, 50), rep(2, 50))
x_dist <- dist(x)
x_distmat <- matrix(0, nrow = 100, ncol = 100)
pnt = 0
for(i in 1:99){
  len = 100-i
  x_distmat[(i+1):100, i] = x_dist[(pnt+1):(pnt+len)]
  pnt = pnt+len
}
x_distmat = x_distmat + t(x_distmat)
x_smds <- list(
  a.0 =  TrainSuperMDS(d = x_distmat, y = y_lab, alpha = 0.0),
  a.2 =  TrainSuperMDS(d = x_distmat, y = y_lab, alpha = 0.2),
  a.5 =  TrainSuperMDS(d = x_distmat, y = y_lab, alpha = 0.5),
  a.10 = TrainSuperMDS(d = x_distmat, y = y_lab, alpha = 1)
)

pdf("sim_smds.pdf", height = 2.5, width = 10)
par(mfrow = c(1,4), mar = c(2,2,2,1))
plot(x_smds$a.0$z, col = y_lab, 
     pch = ifelse(y_lab == 1, 16, 17), cex = 1.5,
     xlab = "", ylab = "", main = "alpha = 0")
plot(x_smds$a.2$z, col = y_lab, 
     pch = ifelse(y_lab == 1, 16, 17), cex = 1.5,
     xlab = "", ylab = "", main = "alpha = 0.2")
plot(x_smds$a.5$z, col = y_lab, 
     pch = ifelse(y_lab == 1, 16, 17), cex = 1.5,
     xlab = "", ylab = "", main = "alpha = 0.5")
plot(x_smds$a.10$z, col = y_lab, 
     pch = ifelse(y_lab == 1, 16, 17), cex = 1.5,
     xlab = "", ylab = "", main = "alpha = 1")
dev.off()