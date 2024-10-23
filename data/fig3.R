library(ggplot2)
source('fig_util.R')

method_v <- c('fmds', 'mds', 'smds', 'umap_s', 'umap_u', 'tsne', 'iso')
params1_v <- (1:4)/4
params2_v <- c(5, 10, 20, 30)

eval_df <- data.frame(matrix(ncol=6, nrow=0))
colnames(eval_df) <- c('method', 'hyperparameter', 'stress', 'corr_dist', 'p_ratio', 'corr_f')

sim_data <- as.matrix(read.csv('result/Evaluation/sim_1-data.csv'))
y_sim <- as.matrix(read.csv('result/Evaluation/sim_1-Y.csv'))

for(m in method_v){
  if(m=='mds'){
    z_embed <- read.csv('result/Evaluation/sim_1-mds-Z.csv')
    res <- rp_eval2(dm=as.matrix(dist(sim_data)), y_orig=y_sim, z_emb=z_embed)
    s <- get_stress(x_dist=dist(sim_data), z_dist=dist(z_embed))
    corr_dist <- get_pearson_corr(x_dist=dist(sim_data), z_dist=dist(z_embed))
    eval_df[nrow(eval_df)+1,] <- list(m, p, s, corr_dist, res$q, res$rho)
    next
  }
  params_v <- if(m=='fmds' | m=='smds') params1_v else params2_v
  for(p in params_v){
    z_embed <- 
      if(m=='fmds' | m=='smds'){
        read.csv(sprintf('result/Evaluation/sim_1-%s-%.2f-Z.csv', m, p))
      } else read.csv(sprintf('result/Evaluation/sim_1-%s-%02d-Z.csv', m, p))
    res <- rp_eval2(dm=as.matrix(dist(sim_data)), y_orig=y_sim, z_emb=z_embed)
    s <- get_stress(x_dist=dist(sim_data), z_dist=dist(z_embed))
    corr_dist <- get_pearson_corr(x_dist=dist(sim_data), z_dist=dist(z_embed))
    eval_df[nrow(eval_df)+1,] <- list(m, p, s, corr_dist, res$q, res$rho)
  }
}


# Panel A: stress v. corr_dist
ggplot(data=eval_df) +
  geom_point(aes(x=stress, y=corr_dist, shape=method, color=method, 
                 fill=as.factor(hyperparameter)), size=2, stroke=0.35) +
  scale_shape_manual(values=c(21,3,16,22,4,24,25)) + 
    # ('fmds', 'iso', 'mds', 'smds', 'tsne', 'umap_s', 'umap_u')
  scale_color_manual(values=c('red','black','black','black','black','black',
                              'black','black')) +
  scale_fill_manual(values=c('grey99', 'grey75', 'grey50', 'grey25', 'grey99', 
                             'grey75', 'grey50', 'grey25')) +
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

ggsave('result/fig3A.pdf', width=2.1, height=2, units='in')


# Panel B: p_ratio v. corr_f
ggplot(data=eval_df) +
  geom_point(aes(x=corr_f, y=p_ratio, shape=method, color=method, 
                 fill=as.factor(hyperparameter)), size=2, stroke=0.35) +
  scale_shape_manual(values=c(21,3,16,22,4,24,25)) + 
    # ('fmds', 'iso', 'mds', 'smds', 'tsne', 'umap_s', 'umap_u')
  scale_color_manual(values=c('red','black','black','black','black','black',
                              'black','black')) +
  scale_fill_manual(values=c('grey99', 'grey75', 'grey50', 'grey25', 'grey99', 
                             'grey75', 'grey50', 'grey25')) +
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

ggsave('result/fig3B.pdf', width=3, height=2, units='in')