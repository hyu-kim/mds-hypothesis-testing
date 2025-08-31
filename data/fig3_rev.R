library(ggplot2)
library(vegan)
library('dplyr')
source('fig_util.R')

# custom plot option
myplot <- function(x, y, ...){
  plot(x, y, xlab = "", ylab = "", xaxt = 'n', yaxt ='n',
       # xaxp = c(0, 0.8, 2),
       xlim = c(0, 1), ylim = c(0, 1.6),
       col = alpha("black", 0.1), 
       pch=16, cex=0.25, cex.axis=1, 
       tcl=-0.2, lwd=0.5, ...
  )
  axis(side=1, at=c(0, 0.5, 1), lwd=0.5, tck=-0.05)
  axis(side=2, lwd=0.5, tck=-0.05, las=2)
}


## import and process
# v_lambda = (0:20)/20
v_lambda1 = c(0:5)/5
v_lambda2 = c(c(0:5)/50, c(4,8,12,16,20)/20)
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

df_eval_stat <- df_eval %>%
  group_by(method, lambda) %>%
  summarize(stress_mean = mean(stress), stress_std = sd(stress),
            p_diff_mean = mean(p_diff), p_diff_std = sd(p_diff),
            pearson_corr_mean = mean(pearson_corr), pearson_corr_std = sd(pearson_corr),
            n = n())


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


## A. Shepard plot
x_dist <- vegdist(t(readRDS('result/ScalingStudy/sim_rev_1/sim_rev_1-N200-data.Rds')), method="bray")
z_dist_fmds_0 <- get_dist_mat(read.csv('result/ScalingStudy/sim_rev_1/sim_rev_1-N200-fmds-0.00-Z.csv'))
z_dist_fmds_1 <- get_dist_mat(read.csv('result/ScalingStudy/sim_rev_1/sim_rev_1-N200-fmds-1.00-Z.csv'))
z_dist_smds_0 <- get_dist_mat(read.csv('result/HyperparameterStudy/SMDS/sim_rev_1-smds-0.00-Z.csv'))
z_dist_smds_1 <- get_dist_mat(read.csv('result/HyperparameterStudy/SMDS/sim_rev_1-smds-1.00-Z.csv'))


pdf("figures/fig3A_rev.pdf", width = 5, height = 1.6)
par(mfrow = c(1, 4), mar = c(1.5,2,0.2,0.3), mgp = c(0,0.5,0), lwd=0.75)
myplot(x_dist, z_dist_fmds_0)
myplot(x_dist, z_dist_fmds_1)
myplot(x_dist, z_dist_smds_0)
myplot(x_dist, z_dist_smds_1)
dev.off()

df_eval_stat <- df_eval_stat[(df_eval_stat$lambda==0)|(df_eval_stat$lambda>=0.2),]

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