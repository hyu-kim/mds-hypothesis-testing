library(ggplot2)
library(vegan)
library('dplyr')
source('fig_util.R')

## A. Shepard plot
x_dist <- vegdist(t(readRDS('result/ScalingStudy/sim_rev_1/sim_rev_1-N200-data.Rds')), method="bray")
z_dist_fmds_0 <- get_dist_mat(read.csv('result/ScalingStudy/sim_rev_1/sim_rev_1-N200-fmds-0.00-Z.csv'))
z_dist_fmds_1 <- get_dist_mat(read.csv('result/ScalingStudy/sim_rev_1/sim_rev_1-N200-fmds-1.00-Z.csv'))
z_dist_smds_0 <- get_dist_mat(read.csv('result/HyperparameterStudy/SMDS/sim_rev_1-smds-0.00-Z.csv'))
z_dist_smds_1 <- get_dist_mat(read.csv('result/HyperparameterStudy/SMDS/sim_rev_1-smds-1.00-Z.csv'))

viz_df <- data.frame(matrix(nrow=0, ncol=5))
colnames(viz_df) <- c('panel', 'lambda', 'method', 'd_x', 'd_z')

for(l in c(0, 1)){
  x_dist <- vegdist(t(readRDS('result/ScalingStudy/sim_rev_1/sim_rev_1-N200-data.Rds')), method="bray")
  z_dist_fmds <- get_dist_mat(read.csv(sprintf('result/ScalingStudy/sim_rev_1/sim_rev_1-N200-fmds-%.2f-Z.csv', l)))
  z_dist_smds <- get_dist_mat(read.csv(sprintf('result/HyperparameterStudy/SMDS/sim_rev_1-smds-%.2f-Z.csv', l)))
  
  viz_df <- rbind(viz_df, data.frame(panel = factor(2*l + 1), lambda = l, 
                                     method = 'fmds', d_x = as.vector(x_dist),
                                     d_z = as.vector(z_dist_fmds)))
  
  viz_df <- rbind(viz_df, data.frame(panel = factor(2*l + 2), lambda = l, 
                                     method = 'smds', d_x = as.vector(x_dist),
                                     d_z = as.vector(z_dist_smds)))
}

ggplot(viz_df) + 
  geom_point(aes(x=d_x, y=d_z), shape=".", alpha=0.1) +
  facet_wrap2(~panel, scales = "free", nrow = 1) +
  scale_y_continuous(limits = c(0, 1.5), breaks=seq(0,1.5,0.5)) +
  scale_x_continuous(limits = c(0, 1), breaks=seq(0,1,0.5)) +
  theme(strip.background = element_rect(fill=NA),
        strip.text = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(fill = "transparent", color = 'black', size=0.5),
        legend.position = "Bottom",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=8, colour='black'),
        # axis.text.y = element_text(size=8, colour='black'),
        axis.text.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth=0.25, colour = 'black')
  )

ggsave('figures/fig3A_rev.pdf', width=4.8, height=1.9, units='in')


## B, C. import and process
v_lambda1 = c(0:10)/10
df_eval <- data.frame(matrix(ncol=6, nrow=0))
colnames(df_eval) <- c('method', 'rep', 'lambda', 'stress', 'p_diff', 'pearson_corr')

for(rep in 1:3){
  x_df <- readRDS(sprintf('result/ScalingStudy/sim_rev_%g/sim_rev_%g-N200-data.Rds', rep, rep))
  x_dist <- vegdist(t(x_df), method="bray")
  y_df <- read.csv(sprintf('result/ScalingStudy/sim_rev-N200-Y.csv'))
  
  for(lambda in v_lambda1){ # FMDS
    z_lambda_df <- read.csv(sprintf('result/ScalingStudy/sim_rev_%g/sim_rev_%g-N200-fmds-%.2f-Z.csv', rep, rep, lambda))
    z_dist <- get_dist_mat(z_lambda_df)
    stress_lambda <- get_stress(x_dist, z_dist)
    # stress_raw_lambda <- get_stress_raw(x_dist, z_dist)
    p_diff <- get_p_diff(z=z_lambda_df, d=x_dist, y=y_df)
    corr_lambda <- get_pearson_corr(x_dist, z_dist)
    df_eval[nrow(df_eval)+1,] <- list('FMDS', rep, lambda, stress_lambda, p_diff, corr_lambda)
  }
  
  for(lambda in v_lambda1){ # SMDS
    z_lambda_df <- read.csv(sprintf('result/HyperparameterStudy/SMDS/sim_rev_%d-smds-%.2f-Z.csv',rep,lambda))
    z_dist <- get_dist_mat(z_lambda_df)
    stress_lambda <- get_stress(x_dist, z_dist)
    # stress_raw_lambda <- get_stress_raw(x_dist, z_dist)
    p_diff <- get_p_diff(z=z_lambda_df, d=x_dist, y=y_df)
    corr_lambda <- get_pearson_corr(x_dist, z_dist)
    df_eval[nrow(df_eval)+1,] <- list('SMDS', rep, lambda, stress_lambda, p_diff, corr_lambda)
  }
}

# save or load objects
saveRDS(df_eval, "result/fig3_rev_dfeval.Rds")
df_eval <- readRDS("result/fig3_rev_dfeval.Rds")

df_eval_stat <- df_eval %>%
  group_by(method, lambda) %>%
  summarize(stress_mean = mean(stress), stress_std = sd(stress),
            p_diff_mean = mean(p_diff), p_diff_std = sd(p_diff),
            pearson_corr_mean = mean(pearson_corr), pearson_corr_std = sd(pearson_corr),
            n = n())

# df_eval_stat <- df_eval_stat[(df_eval_stat$lambda==0)|(df_eval_stat$lambda>=0.2),]

## B. Lambda v dist-corr
ggplot(data=df_eval_stat) +
  geom_point(aes(x=lambda, y=pearson_corr_mean, shape=method), size=2, stroke=0.4, color='black') +
  geom_errorbar(aes(x=lambda, ymin=pearson_corr_mean-pearson_corr_std, 
                    ymax=pearson_corr_mean+pearson_corr_std),
                width=0.05, color='black', size=0.25) +
  scale_shape_manual(values=c(1,2)) +
  # scale_y_continuous(limits=c(0.25,1), breaks = (1:5)/5) +
  scale_x_continuous(breaks = (0:5)/5) +
  theme(strip.background = element_rect(fill=NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(fill = "transparent", color = 'black', size=0.5),
        legend.position = "None",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=8, colour='black'),
        axis.text.y = element_text(size=8, colour='black'),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth=0.5)
  )

ggsave('figures/fig3B_rev.pdf', width=2.0, height=1.95, units='in')


## C. Lambda vs Stress-1
ggplot(data=df_eval_stat) +
  geom_point(aes(x=lambda, y=stress_mean, shape=method), size=2, stroke=0.4, color='black') +
  geom_errorbar(aes(x=lambda, ymin=stress_mean-stress_std, ymax=stress_mean+stress_std),
                width=0.05, color='black', size=0.25) +
  scale_shape_manual(values=c(1,2)) +
  scale_y_continuous(breaks = (0:5)/5) +
  scale_x_continuous(breaks = (0:5)/5) +
  theme(strip.background = element_rect(fill=NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(fill = "transparent", color = 'black', size=0.5),
        legend.position = "None",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=8, colour='black'),
        axis.text.y = element_text(size=8, colour='black'),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth=0.5)
  )

ggsave('figures/fig3C_rev.pdf', width=2.0, height=1.95, units='in')


## Numbers to insert in main text
df_eval_rate <- data.frame(matrix(ncol=3, nrow=0))
colnames(df_eval_rate) <- c('method', 'rep', 'stress_rate')
for(m in c('FMDS', 'SMDS')){
  for(r in 1:3){
    df_eval_rate[nrow(df_eval_rate)+1,] <- 
      list(m, r, 
           df_eval$stress[df_eval$method==m & df_eval$rep==r & df_eval$lambda==1] -
             df_eval$stress[df_eval$method==m & df_eval$rep==r & df_eval$lambda==0]
      )
  }
}

print(mean(df_eval_rate$stress_rate[df_eval_rate$method=='SMDS']) / 
        mean(df_eval_rate$stress_rate[df_eval_rate$method=='FMDS']))