library(tmvtnorm)
source('mm.R')
N <- 100
del <- 1/(2*sqrt(2))*0.1

n1 <- as.matrix(c(1,1,1,1))/2
n2 <- as.matrix(c(1,1,1,-3))/sqrt(12)
n3 <- as.matrix(c(1,1,-2,0))/sqrt(6)
n4 <- as.matrix(c(1,-1,0,0))/sqrt(2)
lambda_diag <- matrix(c(0.0001,0,0,0, 0,0.04,0,0, 0,0,0.04,0, 0,0,0,0.01), nrow=4)

covmat <- round(cbind(n1,n2,n3,n4) %*% lambda_diag %*% t(cbind(n1,n2,n3,n4)),10)

sim_data <- list(
  data = rbind(
    rtmvnorm(N/2, c(0.25+del,0.25-del,0.25,0.25), sigma = covmat, lower=rep(0,4), upper=rep(1,4)),
    rtmvnorm(N/2, c(0.25-del,0.25+del,0.25,0.25), sigma = covmat, lower=rep(0,4), upper=rep(1,4))
    )
)

sim_data$data <- apply(sim_data$data, 1, function(i) i/sum(i))
sim_data$data <- t(sim_data$data)

sim_data$data <- data.frame(sim_data$data)
sim_data$Y <- c(rep(1, N/2), rep(2, N/2))
sim_data$distmat <- as.matrix(dist(sim_data$data[, 1:4]))
z0 <- cmdscale(dist(sim_data$data[,1:4]), k = 2)

write.csv(sim_data$data, 'result/HyperparameterStudy/sim-data.csv', row.names=FALSE)
write.csv(sim_data$Y, 'result/HyperparameterStudy/sim-Y.csv', row.names=FALSE)