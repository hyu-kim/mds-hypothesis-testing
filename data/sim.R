#3. Simulation
library(MASS)
library(parallel)
library(pbmcapply)
library(vegan)
source('mm.R')

## Duplicates and hyperparameters
# single run of MM
run_mm_cmds <- function(params) {
  r <- params$replicate
  N <- params$size
  l <- params$lambda
  
  data_mat <- readRDS(sprintf('result/ScalingStudy/sim_rev2_%d/sim_rev2_%d-N%d-data.Rds', r, r, N))
  dist_mat <- as.matrix(vegdist(t(data_mat), method="bray"))
  y_df <- as.matrix(read.csv(sprintf('result/ScalingStudy/sim_rev2-N%d-Y.csv', N)))
  
  print(sprintf("Data Read. N: %d, Replicate: %d, For lambda: %g", N, r, l))
  
  z0 <- cmdscale(dist_mat, k = 2)
  res <- mm_cmds(nit = 50, lambda = l, z0 = z0, D = dist_mat, y = y_df,
                 threshold = 0.1, dataset = sprintf('sim_rev2_%d_N%d', r, N),
                 folder = 'result/Revision2_Dec2025')
  
  print(sprintf("MM Done. N: %d, Replicate: %d, For lambda: %g", N, r, l))
}

# dataframe of combinations
replicate_list <- c(1)
size_list <- c(500, 500, 500)
lambda_list <- c(0.04, 0.2, 0.8)

param_combinations <- data.frame(replicate = replicate_list, size = size_list, lambda = lambda_list)
results <- pbmclapply(split(param_combinations, seq_len(nrow(param_combinations))),
                      run_mm_cmds,
                      mc.cores = detectCores() - 1)


# 50 iteration
for(N in c(50,100,200,500)){
  print(sprintf("Data Size: %d",N))
  data_mat <- readRDS(sprintf('result/ScalingStudy/sim_rev_3/sim_rev_3-N%d-data.Rds', N))
  # data_mat <- readRDS(sprintf('result/ScalingStudy/sim_rev_%d-N200-data.Rds',r))
  dist_mat <- as.matrix(vegdist(t(data_mat), method="bray"))
  # y_df <- as.matrix(read.csv('result/ScalingStudy/sim_rev-N200-Y.csv'))
  y_df <- as.matrix(read.csv(sprintf('result/ScalingStudy/sim_rev-N%d-Y.csv', N)))
  # dist_mat <- as.matrix(dist(data_df[,1:4]))
  z0 <- cmdscale(dist_mat, k = 2)
  res <- pbmclapply(1:19, function(i){
    x <- c(2:20)[i] * 0.05
    return(mm_cmds(nit = 50, lambda = x, z0 = z0, D = dist_mat, y = y_df,
                   dataset = paste('sim_rev_3-N',N,sep='')))
  }, mc.cores = 128)
  write.csv(data.frame(n_iter=50, dataset=paste('sim_rev_1-N',N,sep='')),
            sprintf('result/ScalingStudy/sim_rev_1_%d-config.csv',N), row.names=FALSE)
  print(sprintf("Size %d done",N))
}


# 200 iteration
for(r in 1:3){
  print(sprintf("Set: %d",r))
  data_df <- as.matrix(read.csv(sprintf('result/HyperparameterStudy/sim_%d-data.csv',r)))
  y_df <- as.matrix(read.csv(sprintf('result/HyperparameterStudy/sim_%d-Y.csv',r)))
  dist_mat <- as.matrix(dist(data_df[,1:4]))
  z0 <- cmdscale(dist_mat, k = 2)
  res <- pbmclapply(1:3, function(i){
    x <- c(0:2)[i] * 0.05
    return(mm_cmds(nit = 20, lambda = x, z0 = z0, D = dist_mat, y = y_df,
                   dataset = paste('sim',r,sep='_')))
  }, mc.cores = 128)
  # res <- mm_cmds(nit = 200, lambda = 0.10, z0 = z0, D = dist_mat, y = y_df, dataset = paste('sim',r,sep='_'))
  print(sprintf("Set %d done",r))
}