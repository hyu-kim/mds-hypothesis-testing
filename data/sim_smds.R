library(parallel)
library(pbmcapply)
library(superMDS)
source('mm.R')

for(r in 1:3){
  print(sprintf("Set: %d",r))
  data_df <- as.matrix(read.csv(sprintf('result/HyperparameterStudy/sim_%d-data.csv',r)))
  y_df <- as.matrix(read.csv(sprintf('result/HyperparameterStudy/sim_%d-Y.csv',r)))
  dist_mat <- as.matrix(dist(data_df[,1:4]))
  z0 <- cmdscale(dist_mat, k = 2)
  pbmclapply(1:6, function(i){
    x <- c(1:4)[i] * 0.02
    res <- TrainSuperMDS(d = dist_mat, y = y_df, alpha = x)
    write.csv(res$z, sprintf('result/HyperparameterStudy/SMDS/sim_%d-smds-%.2f-Z.csv', 
                             r, x), row.names=FALSE)
  }, mc.cores = 128)
  print(sprintf("Set %d done",r))
}