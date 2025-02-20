library(ggplot2)
library('dplyr')
source('fig_util.R')

# custom plot option
myplot <- function(x, y, ...){
  plot(x, y, xlab = "", ylab = "", xaxt = 'n', yaxt ='n',
       col = alpha("black", 0.15), 
       xlim = c(0, 0.9), ylim = c(0, 0.9),
       pch=16, cex=0.3, cex.axis=1, 
       tcl=-0.2, lwd=0.5, ...
  )
  axis(side=1, lwd=0.5, tck=-0.05)
  axis(side=2, lwd=0.5, tck=-0.05, las=2)
}


## import and process
# v_lambda = (0:20)/20
v_lambda1 = c(0,1,2,3,4,8,12,16,20)/20
v_lambda2 = c(c(0:5)/50, c(4,8,12,16,20)/20)
df_eval <- data.frame(matrix(ncol=6, nrow=0))
colnames(df_eval) <- c('method', 'rep', 'lambda', 'stress', 'p_diff', 'pearson_corr')

for(rep in 1:3){
  for(lambda in v_lambda1){ # FMDS
    # lambda <- v_lambda[i]
    if(lambda==0.15) {
      z_lambda_df <- read.csv(sprintf('result/HyperparameterStudy/iter_200/sim_%d-fmds-%.2f-Z.csv', rep, lambda))
      } else z_lambda_df <- read.csv(sprintf('result/HyperparameterStudy/sim_%d/sim_%d-fmds-%.2f-Z.csv', rep, rep, lambda))
    y_df <- read.csv(sprintf('result/HyperparameterStudy/sim_%d-Y.csv', rep))
    x_df <- read.csv(sprintf('result/HyperparameterStudy/sim_%d-data.csv', rep))
    x_dist <- get_dist_mat(x_df)
    z_dist <- get_dist_mat(z_lambda_df)
    stress_lambda <- get_stress(x_dist, z_dist)
    # stress_raw_lambda <- get_stress_raw(x_dist, z_dist)
    p_diff <- get_p_diff(z=z_lambda_df, d=x_dist, y=y_df)
    corr_lambda <- get_pearson_corr(x_dist, z_dist)
    df_eval[nrow(df_eval)+1,] <- list('FMDS', rep, lambda, stress_lambda, p_diff, corr_lambda)
  }
  
  for(lambda in v_lambda2){ # SMDS
    z_lambda_df <- read.csv(sprintf('result/HyperparameterStudy/SMDS/sim_%d-smds-%.2f-Z.csv',rep,lambda))
    x_dist <- get_dist_mat(x_df)
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


## A. Lambda v p_diff
ggplot(data=df_eval_stat) +
  geom_point(aes(x=lambda, y=p_diff_mean, shape=method), size=1.5, stroke=0.25) +
  geom_errorbar(aes(x=lambda, ymin=p_diff_mean-p_diff_std, ymax=p_diff_mean+p_diff_std),
                width=0.02, color='black', size=0.25) +
  scale_shape_manual(values=c(1,2)) +
  # scale_y_continuous(limits=c(0.25,1), breaks = (1:5)/5) +
  scale_x_continuous(trans = 'log10') +
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

ggsave('figures/fig2A_v2.pdf', width=1.6, height=1.55, units='in')

## A. Lambda v dist-corr
ggplot(data=df_eval_stat) +
  geom_point(aes(x=lambda, y=pearson_corr_mean, shape=method), size=1.5, stroke=0.25) +
  geom_errorbar(aes(x=lambda, ymin=pearson_corr_mean-pearson_corr_std, 
                    ymax=pearson_corr_mean+pearson_corr_std),
                width=0.05, color='black', size=0.25) +
  scale_shape_manual(values=c(1,2)) +
  scale_y_continuous(limits=c(0.25,1), breaks = (1:5)/5) +
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

ggsave('figures/fig2A.pdf', width=1.6, height=1.55, units='in')


## B. Lambda vs Stress-1
ggplot(data=df_eval_stat) +
  geom_point(aes(x=lambda, y=stress_mean, shape=method), size=1.5, stroke=0.25) +
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

ggsave('figures/fig2B.pdf', width=1.6, height=1.55, units='in')


## C. Lambda v raw stress
stress_raw_max <- max(df_eval$stress_raw)
ggplot(data=df_eval_stat) +
  geom_point(aes(x=lambda, y=stress_raw_mean/1, shape=method), size=1.5, stroke=0.25) +
  geom_errorbar(aes(x=lambda, ymin=(stress_raw_mean-stress_raw_std)/1, 
                    ymax=(stress_raw_mean+stress_raw_std)/1),
                width=0.05, color='black', size=0.25) +
  scale_shape_manual(values=c(1,2)) +
  # scale_y_continuous(limits=c(0,1), breaks = (0:5)/5) +
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
        axis.ticks = element_line(linewidth=0.25, colour = 'black')
  )

ggsave('figures/fig2C.pdf', width=1.6, height=1.55, units='in')


## D. Shepard plot
x_dist <- get_dist_mat(read.csv('result/HyperparameterStudy/sim_1-data.csv'))
z_dist_fmds_p0 <- get_dist_mat(read.csv('result/HyperparameterStudy/sim_1/sim_1-fmds-0.00-Z.csv'))
z_dist_fmds_p6 <- get_dist_mat(read.csv('result/HyperparameterStudy/sim_1/sim_1-fmds-0.60-Z.csv'))
z_dist_smds_p0 <- get_dist_mat(read.csv('result/HyperparameterStudy/SMDS/sim_1-smds-0.00-Z.csv'))
z_dist_smds_p6 <- get_dist_mat(read.csv('result/HyperparameterStudy/SMDS/sim_1-smds-0.60-Z.csv'))


pdf("figures/fig2D.pdf", width = 5.1, height = 1.1)
par(mfrow = c(1, 4), mar = c(1.5,3,0.2,0.2), mgp = c(0,0.6,0), lwd=0.75)
myplot(x_dist, z_dist_fmds_p0)
myplot(x_dist, z_dist_smds_p0)
myplot(x_dist, z_dist_fmds_p6)
myplot(x_dist, z_dist_smds_p6)
dev.off()