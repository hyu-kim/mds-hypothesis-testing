## computes F-MDS w/o convergence limit
library(vegan)
library(parallel)
source("mm.R")

run_mm_cmds <- function(params) {
  require(vegan)
  require("mm.R")
  require("permanova_with_config.R")
  r <- params$replicate
  N <- params$size
  l <- params$lambda
  
  data_mat <- readRDS('result/ScalingStudy/sim_rev_1/sim_rev_1-N100-data.Rds')
  dist_mat <- as.matrix(vegdist(t(data_mat), method="bray"))
  y_df <- as.matrix(read.csv(sprintf('result/ScalingStudy/sim_rev-N100-Y.csv')))
  
  print(sprintf("Data Read. N: %d, Replicate: %d, For lambda: %g", N, r, l))
  
  z0 <- cmdscale(dist_mat, k = 2)
  res <- mm_cmds(nit = 50, lambda = l, z0 = z0, D = dist_mat, y = y_df,
                 threshold = 0, dataset = 'sim_rev_1-N100', 
                 folder = "result/Revision2_Dec2025")
  
  write.csv(res, sprintf('Simulated/F-MDS/Results/sim_rev_%d-N%d-fmds-%.2f-Z', r, N, l), row.names=FALSE)
  
  print(sprintf("MM Done. N: %d, Replicate: %d, For lambda: %g", N, r, l))
}

# dataframe of combinations
replicate_list <- c(1)
size_list <- c(100)
lambda_list <- c(0.2, 0.5, 0.8)

# loop and run parellel
param_combinations <- expand.grid(replicate_list, size_list, lambda_list)
colnames(param_combinations) <- c("replicate", "size", "lambda")

num_cores <- detectCores()
cl <- makeCluster(num_cores - 1) # Leave one core free for system tasks

clusterExport(cl, ls(envir = .GlobalEnv), envir = environment())

results <- parLapply(cl, split(param_combinations, seq_len(nrow(param_combinations))), run_mm_cmds)

stopCluster(cl)
