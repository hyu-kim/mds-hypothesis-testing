# configuration by pre-designed basis
n1 <- as.matrix(c(1,1,1,1))/2
n2 <- as.matrix(c(1,1,1,-3))/sqrt(12)
n3 <- as.matrix(c(1,1,-2,0))/sqrt(6)
n4 <- as.matrix(c(1,-1,0,0))/sqrt(2)
n_mat <- cbind(n1, n2, n3, n4)

sim_data_mat <- as.matrix(read.csv(sprintf('result/HyperparameterStudy/sim_%d-data.csv',1)))
sim_data_norm <- apply(sim_data_mat, 1, function(i) i/sum(i))
sim_data_norm <- t(sim_data_norm)
y <- as.matrix(read.csv(sprintf('result/HyperparameterStudy/sim_%d-Y.csv',1)))

sim_proj_mat <- sim_data_norm %*% n_mat

pdf("figures/Fig_S1.pdf", width = 6.5, height = 6.5)

par(mfrow = c(3, 3), mar = c(3,3,1,1), mgp = c(2,1,0))

for(r in 1:3){
  plot(sim_proj_mat[,1], sim_proj_mat[,4],
       xlab = "v1", ylab = "v4", 
       ylim = c(-0.4, 0.4),
       pch = c(16,17)[y], col = c('black', 'red')[y])
  plot(sim_proj_mat[,2], sim_proj_mat[,4],
       xlab = "v2", ylab = "v4", 
       xlim = c(-0.4, 0.4), ylim = c(-0.4, 0.4),
       pch = c(16,17)[y], col = c('black', 'red')[y])
  plot(sim_proj_mat[,3], sim_proj_mat[,4],
       xlab = "v3", ylab = "v4", 
       xlim = c(-0.4, 0.4), ylim = c(-0.4, 0.4),
       pch = c(16,17)[y], col = c('black', 'red')[y])
  
  plot(sim_proj_mat[,1], sim_proj_mat[,3],
       xlab = "v1", ylab = "v3",
       ylim = c(-0.4, 0.4),
       pch = c(16,17)[y], col = c('black', 'red')[y])
  plot(sim_proj_mat[,2], sim_proj_mat[,3],
       xlab = "v2", ylab = "v3",
       xlim = c(-0.4, 0.4), ylim = c(-0.4, 0.4),
       pch = c(16,17)[y], col = c('black', 'red')[y])
  plot.new()
  
  plot(sim_proj_mat[,1], sim_proj_mat[,2],
       xlab = "v1", ylab = "v2",
       ylim = c(-0.4, 0.4),
       pch = c(16,17)[y], col = c('black', 'red')[y])
  plot.new()
  plot.new()
}

dev.off()