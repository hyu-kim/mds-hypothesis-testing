library(ggplot2)
library(vegan)
library('dplyr')
source('fig_util.R')


## import and process
v_lambda1 = c(0:10)/10
df_eval <- data.frame(matrix(ncol=7, nrow=0))
colnames(df_eval) <- c('N', 'method', 'rep', 'lambda', 'stress', 'p_diff', 'pearson_corr')

for(N in c(50, 100, 500)){
  for(rep in 1:3){
    message("N: ", N, ", rep: ", rep)
    x_df <- readRDS(sprintf('result/ScalingStudy/sim_rev_%g/sim_rev_%g-N%g-data.Rds', rep, rep, N))
    x_dist <- vegdist(t(x_df), method="bray")
    y_df <- read.csv(sprintf('result/ScalingStudy/sim_rev-N%g-Y.csv', N))
    
    for(lambda in v_lambda1){ # FMDS
      z_lambda_df <- read.csv(sprintf('result/ScalingStudy/sim_rev_%g/sim_rev_%g-N%g-fmds-%.2f-Z.csv',
                                      rep, rep, N, lambda))
      z_dist <- get_dist_mat(z_lambda_df)
      stress_lambda <- get_stress(x_dist, z_dist)
      # stress_raw_lambda <- get_stress_raw(x_dist, z_dist)
      p_diff <- get_p_diff(z=z_lambda_df, d=x_dist, y=y_df)
      corr_lambda <- get_pearson_corr(x_dist, z_dist)
      df_eval[nrow(df_eval)+1,] <- list(N, 'FMDS', rep, lambda, stress_lambda, p_diff, corr_lambda)
    }
    
    for(lambda in v_lambda1){ # SMDS
      z_lambda_df <- read.csv(sprintf('result/HyperparameterStudy/SMDS/sim_rev_%d-N%g-smds-%.2f-Z.csv',
                                      rep, N, lambda))
      z_dist <- get_dist_mat(z_lambda_df)
      stress_lambda <- get_stress(x_dist, z_dist)
      # stress_raw_lambda <- get_stress_raw(x_dist, z_dist)
      p_diff <- get_p_diff(z=z_lambda_df, d=x_dist, y=y_df)
      corr_lambda <- get_pearson_corr(x_dist, z_dist)
      df_eval[nrow(df_eval)+1,] <- list(N, 'SMDS', rep, lambda, stress_lambda, p_diff, corr_lambda)
    }
  }
}

# save or load objects
saveRDS(df_eval, "result/S3Fig_rev_dfeval.Rds")
df_eval <- readRDS("result/S3Fig_rev_dfeval.Rds")

df_eval_stat <- df_eval %>%
  group_by(N, method, lambda) %>%
  summarize(stress_mean = mean(stress), stress_std = sd(stress),
            p_diff_mean = mean(p_diff), p_diff_std = sd(p_diff),
            pearson_corr_mean = mean(pearson_corr), pearson_corr_std = sd(pearson_corr),
            n = n())

# df_eval_stat <- df_eval_stat[(df_eval_stat$lambda==0)|(df_eval_stat$lambda>=0.2),]


## A. Lambda v dist-corr
labels_v <- c("50" = 'N = 50', "100" = 'N = 100', "500" = 'N = 500')

ggplot(data=df_eval_stat) +
  geom_point(aes(x=lambda, y=pearson_corr_mean, shape=method), size=2, stroke=0.4, color='black') +
  geom_errorbar(aes(x=lambda, ymin=pearson_corr_mean-pearson_corr_std, 
                    ymax=pearson_corr_mean+pearson_corr_std),
                width=0.05, color='black', size=0.25) +
  scale_shape_manual(values=c(1,2), labels = c("F-MDS", "SuperMDS")) +
  scale_x_continuous(breaks = (0:5)/5) +
  facet_wrap(~N, axes="all", nrow=1, labeller = as_labeller(labels_v)) +
  labs(x = "Hyperparameter", y = "Pearson correlation") +
  theme(strip.background = element_rect(fill=NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(fill = "transparent", color = 'black', size=0.5),
        legend.position = "bottom",
        # axis.title.x = element_blank(),
        # axis.title.y = element_blank(),
        # axis.text.x = element_text(size=9, colour='black'),
        # axis.text.y = element_text(size=9, colour='black'),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth=0.5)
  )

ggsave('figures/Fig_S3A_rev.pdf', width=6.5, height=3, units='in')


## B. Lambda vs Stress-1
ggplot(data=df_eval_stat) +
  geom_point(aes(x=lambda, y=stress_mean, shape=method), size=2, stroke=0.4, color='black') +
  geom_errorbar(aes(x=lambda, ymin=stress_mean-stress_std, ymax=stress_mean+stress_std),
                width=0.05, color='black', size=0.25) +
  scale_shape_manual(values=c(1,2), labels = c("F-MDS", "SuperMDS")) +
  scale_y_continuous(breaks = (0:5)/5) +
  scale_x_continuous(breaks = (0:5)/5) +
  facet_wrap(~N, axes="all", nrow=1, labeller = as_labeller(labels_v)) +
  labs(x = "Hyperparameter", y = "Stress-1") +
  theme(strip.background = element_rect(fill=NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(fill = "transparent", color = 'black', size=0.5),
        legend.position = "bottom",
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth=0.5)
  )

ggsave('figures/Fig_S3B_rev.pdf', width=6.5, height=3, units='in')