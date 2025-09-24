### Computes grid search objective function for hyperparameter selection
library(vegan)
library(ggplot2)
library(ggforce)
library(cowplot)
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

df_eval$'f_obj' <- log(df_eval$n_epoch) * (1-df_eval$pearson_corr)

saveRDS(df_eval, "result/S3Appendix_dfeval.Rds")
df_eval <- readRDS("result/S3Appendix_dfeval.Rds")


# find arg minimum lambda
argmin_df_eval <- data.frame(matrix(ncol=4, nrow=0))
colnames(argmin_df_eval) <- c('N', 'rep', 'lambda', 'f_obj')

for(N in c(50, 100, 200, 500)){
  for(rep in 1:3){
    df_eval_f <- df_eval[df_eval$N==N & df_eval$rep==rep, ]
    lambda <- df_eval_f$lambda[which.min(df_eval_f$f_obj)]
    f_obj <- min(df_eval_f$f_obj)
    argmin_df_eval[nrow(argmin_df_eval)+1,] <- list(N, rep, lambda, f_obj)
  }
}


# statistical test
result <- cor.test(argmin_df_eval$N, argmin_df_eval$lambda) # spearman
print(result)
print(mean(argmin_df_eval$lambda))
print(sd(argmin_df_eval$lambda))


labels_v <- c("50" = 'N = 50', "100" = 'N = 100', "200" = 'N = 200', "500" = 'N = 500')

p1 <- 
  ggplot() +
    geom_line(data = df_eval, aes(x=lambda, y=f_obj, color=as.factor(rep)), linewidth=0.25) +
    geom_point(data = df_eval, aes(x=lambda, y=f_obj, shape=as.factor(rep)), size=1.5, stroke=0.5, color='black') +
    geom_point(data = argmin_df_eval, aes(x=lambda, y=f_obj), shape=4, color='red', size=2, stroke=0.75) +
    facet_wrap(~N, axes="all", nrow=1, labeller = as_labeller(labels_v), scales = "free") +
    scale_shape_manual(values=c(0,1,2)) + 
    scale_color_manual(values=c('black', 'black', 'black')) +
    labs(x = expression(lambda), y = bquote(f[obj] (lambda))) +
    theme(strip.background = element_rect(fill=NA),
          panel.background = element_rect(fill = "transparent", color = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank(),
          panel.border = element_rect(fill = "transparent", color = 'black', size=0.5),
          legend.position = "none",
          legend.title=element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_line(linewidth=0.25, colour = 'black'),
          text = element_text(size = 11)
    )

p2 <- 
  ggplot() +
    geom_sina(data = argmin_df_eval, aes(x=as.factor(N), y=lambda, shape=as.factor(rep)), color='black', size=2.5, stroke=0.5, maxwidth = 0.5) +
    facet_wrap(~N, axes="all", nrow=1, labeller = as_labeller(labels_v), scales = "free_x") +
    scale_shape_manual(values=c(0,1,2)) +
    scale_y_continuous(limits = c(0,1)) +
    scale_x_discrete(labels= labels_v) +
    labs(x = NULL, y = bquote(lambda[min])) +
    theme(strip.background = element_rect(fill=NA),
          strip.text.x = element_blank(),
          panel.background = element_rect(fill = "transparent", color = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank(),
          panel.border = element_rect(fill = "transparent", color = 'black', size=0.5),
          legend.position = "none",
          legend.title=element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_line(linewidth=0.25, colour = 'black'),
          # axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          text = element_text(size = 11)
    )

plot_grid(p1, p2, labels='AUTO', ncol=1, rel_heights = c(1.5, 1))

ggsave('figures/Fig_S4_rev_gridsearch.pdf', width=6.5, height=4, units='in')