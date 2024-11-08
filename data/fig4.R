library(ggplot2)
source('fig_util.R')

method_v <- c('fmds','mds', 'smds', 'umap_s', 'umap_u', 'tsne', 'iso', 'nn')
params1_v <- (1:4)/5
params2_v <- c(5, 10, 20, 30)
params3_v <- c(5, 7, 10)

eval_df <- data.frame(matrix(ncol=6, nrow=0))
colnames(eval_df) <- c('method', 'hyperparameter', 'stress', 'corr_dist', 
                       'p_ratio', 'corr_f')

alga_2_dist <- readRDS('result/Evaluation/alga_2-dist.rds')
alga_2_y <- as.matrix(read.csv('result/Evaluation/alga_2-Y.csv'))

for(m in method_v){
  if(m=='mds'){
    z_embed <- read.csv('result/Evaluation/alga_2-mds-Z.csv')
    res <- rp_eval2(dm=as.matrix(alga_2_dist), y_orig=alga_2_y, z_emb=z_embed)
    # rp_eval2(500, distmat2, y2, nn2[,2:33], nn2$Y)
    s <- get_stress(x_dist=alga_2_dist, z_dist=dist(z_embed))
    corr_dist <- get_pearson_corr(x_dist=alga_2_dist, z_dist=dist(z_embed))
    eval_df[nrow(eval_df)+1,] <- list(m, 0, s, corr_dist, res$q, res$rho)
    next
  }
  else if(m=='nn'){
    alga_2_nn <- readRDS('result/Evaluation/alga_2-nn-data.rds')
    z_embed <- alga_2_nn[,2:33]
    res <- rp_eval2(dm=as.matrix(alga_2_dist), y_orig=alga_2_y, z_emb=z_embed, 
                    y_emb=alga_2_nn$Y)
    s <- get_stress(x_dist=alga_2_dist, z_dist=dist(z_embed))
    corr_dist <- get_pearson_corr(x_dist=alga_2_dist, z_dist=dist(z_embed))
    eval_df[nrow(eval_df)+1,] <- list(m, 0, s, corr_dist, res$q, res$rho)
    next
  }
  params_v <- if(m=='smds' | m=='fmds') params1_v else if (m=='umap_s' | m=='umap_u') params2_v else params3_v
  for(p in params_v){
    z_embed <- 
      if(m=='fmds' | m=='smds'){
        read.csv(sprintf('result/Evaluation/alga_2-%s-%.2f-Z.csv', m, p))
      } else read.csv(sprintf('result/Evaluation/alga_2-%s-%02d-Z.csv', m, p))
    res <- rp_eval2(dm=alga_2_dist, y_orig=alga_2_y, z_emb=z_embed)
    s <- get_stress(x_dist=alga_2_dist, z_dist=dist(z_embed))
    corr_dist <- get_pearson_corr(x_dist=alga_2_dist, z_dist=dist(z_embed))
    eval_df[nrow(eval_df)+1,] <- list(m, p, s, corr_dist, res$q, res$rho)
  }
}


# Panel A: stress v. corr_dist
ggplot(data=eval_df) +
  geom_point(aes(x=stress, y=corr_dist, shape=method, color=method, 
                 fill=as.factor(hyperparameter)), size=2, stroke=0.35) +
  scale_shape_manual(values=c(21,3,16,8,22,4,24,25)) + 
    # ('fmds', 'iso', 'mds', 'nn', 'smds', 'tsne', 'umap_s', 'umap_u')
  scale_color_manual(values=c('red','black', 'black','black','black','black','black','black','black')) +
  scale_fill_manual(values=c('black','grey99', 'grey75', 'grey50', 'grey25', 
                             'grey99', 'grey80', 'grey60', 'grey40', 'grey20')) +
  # scale_y_continuous(breaks = (0:5)/5) +
  # scale_x_continuous(breaks = (0:5)/5) +
  theme(strip.background = element_rect(fill=NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(fill = "transparent", color = 'black', size=0.5),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=8, colour='black'),
        axis.text.y = element_text(size=8, colour='black'),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth=0.25, colour = 'black')
  )

ggsave('result/fig4A.pdf', width=2.1, height=2, units='in')


# Panel B: p_ratio v. corr_f
ggplot(data=eval_df) +
  geom_point(aes(x=corr_f, y=p_ratio, shape=method, color=method, 
                 fill=as.factor(hyperparameter)), size=2, stroke=0.35) +
  scale_shape_manual(values=c(21,3,16,8,22,4,24,25)) + 
    # ('fmds', 'iso', 'mds', 'nn', 'smds', 'tsne', 'umap_s', 'umap_u')
  scale_color_manual(values=c('red','black', 'black','black','black','black','black','black','black')) +
  scale_fill_manual(values=c('black','grey99', 'grey75', 'grey50', 'grey25', 
                             'grey99', 'grey80', 'grey60', 'grey40', 'grey20')) +
  # scale_y_continuous(breaks = (0:5)/5) +
  # scale_x_continuous(breaks = (0:5)/5) +
  scale_y_continuous(breaks = (0:4)/4) +
  theme(strip.background = element_rect(fill=NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(fill = "transparent", color = 'black', size=0.5),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=8, colour='black'),
        axis.text.y = element_text(size=8, colour='black'),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth=0.25, colour = 'black')
  )

ggsave('result/fig4B.pdf', width=3, height=2, units='in')


# Panes: lambda v distance correlation (FMDS, SMDS)
ggplot(data=eval_df[eval_df$method=='fmds'|eval_df$method=='smds',]) +
  geom_point(aes(x=hyperparameter, y=corr_dist, shape=method), size=1.5, stroke=0.25) +
  # geom_errorbar(aes(x=lambda, ymin=pearson_corr_mean-pearson_corr_std, 
  #                   ymax=pearson_corr_mean+pearson_corr_std),
  #               width=0.05, color='black', size=0.25) +
  scale_shape_manual(values=c(1,2)) +
  # scale_y_continuous(limits=c(0.25,1), breaks = (1:5)/5) +
  # scale_x_continuous(breaks = (0:5)/5) +
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