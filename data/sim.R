#3. Simulation
library(MASS)
library(parallel)
library(pbmcapply)
source('mm.R')

# 50 iteration
for(r in 1:3){
  print(sprintf("Set: %d",r))
  data_df <- as.matrix(read.csv(sprintf('result/HyperparameterStudy/sim_%d-data.csv',r)))
  y_df <- as.matrix(read.csv(sprintf('result/HyperparameterStudy/sim_%d-Y.csv',r)))
  dist_mat <- as.matrix(dist(data_df[,1:4]))
  z0 <- cmdscale(dist_mat, k = 2)
  res <- pbmclapply(1:21, function(i){
    x <- c(0:20)[i] * 0.05
    return(mm_cmds(nit = 50, lambda = x, z0 = z0, D = dist_mat, y = y_df,
                   dataset = paste('sim',r,sep='_')))
  }, mc.cores = 128)
  write.csv(data.frame(n_iter=50, dataset=paste('sim',r,sep='_')),
            sprintf('result/HyperparameterStudy/sim_%d/sim_%d-config.csv', r,r), row.names=FALSE)
  print(sprintf("Set %d done",r))
}



# 200 iteration
for(r in 1:3){
  print(sprintf("Set: %d",r))
  data_df <- as.matrix(read.csv(sprintf('result/HyperparameterStudy/sim_%d-data.csv',r)))
  y_df <- as.matrix(read.csv(sprintf('result/HyperparameterStudy/sim_%d-Y.csv',r)))
  dist_mat <- as.matrix(dist(data_df[,1:4]))
  z0 <- cmdscale(dist_mat, k = 2)
  res <- mm_cmds(nit = 200, lambda = 0.15, z0 = z0, D = dist_mat, y = y_df, dataset = paste('sim',r,sep='_'))
  print(sprintf("Set %d done",r))
}