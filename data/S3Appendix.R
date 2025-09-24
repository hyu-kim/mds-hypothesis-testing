### Computes grid search objective function for hyperparameter selection
library(vegan)
library(ggplot2)
source('fig_util.R')

# Load data
df_eval <- data.frame(matrix(ncol=6, nrow=0))
colnames(df_eval) <- c('N', 'rep', 'lambda', 'stress', 'pearson_corr', 'n_epoch')

for(N in c(50, 100, 200, 500)){
  for(rep in 1:3){
    message("N: ", N, ", rep: ", rep)
    x_df <- readRDS(sprintf('result/ScalingStudy/sim_rev_%g/sim_rev_%g-N%g-data.Rds', rep, rep, N))
    x_dist <- vegdist(t(x_df), method="bray")
    y_df <- read.csv(sprintf('result/ScalingStudy/sim_rev-N%g-Y.csv', N))
    
    for(lambda in c(1:10)/10){ # FMDS
      z_lambda_df <- read.csv(sprintf('result/ScalingStudy/sim_rev_%g/sim_rev_%g-N%g-fmds-%.2f-Z.csv',
                                      rep, rep, N, lambda))
      z_dist <- get_dist_mat(z_lambda_df)
      stress_lambda <- get_stress(x_dist, z_dist)
      corr_lambda <- get_pearson_corr(x_dist, z_dist)
      n_iter <- nrow(read.csv(sprintf('result/ScalingStudy/sim_rev_%g/sim_rev_%g-N%d-fmds-%.2f-log.csv',
                                      rep, rep, N, lambda)))
      df_eval[nrow(df_eval)+1,] <- list(N, rep, lambda, stress_lambda, corr_lambda, n_iter)
    }
  }
}

df_eval$'f_obj' <- df_eval$n_epoch * (1-df_eval$pearson_corr)

saveRDS(df_eval, "result/S3Appendix_dfeval.Rds")