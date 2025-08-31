library(vegan)
library(ggplot2)
source('fig_util.R')
source('../source/dreval/continuity.R')
source('../source/dreval/trustworthiness.R')

method_v <- c('fmds', 'mds', 'umap_s', 'umap_u', 'tsne', 'iso')
params1_v <- c(0.2, 0.4, 0.8)
params2_v <- c(5, 10, 50)

eval_df <- data.frame(matrix(ncol=9, nrow=0))
colnames(eval_df) <- c('method', 'hyperparameter', 'k', 'continuity', 
                       'trustworthiness', 'stress', 'corr_dist', 'p_ratio', 'corr_f')

sim_data <- t(readRDS('result/ScalingStudy/sim_rev_1/sim_rev_1-N200-data.Rds'))
rownames(sim_data) <- c(1:200)
sim_data_dist <- vegdist(sim_data, method="bray")
sim_y_mat <- as.matrix(read.csv(sprintf('result/ScalingStudy/sim_rev-N200-Y.csv')))

for(k in c(14,150)){ # local, global
  for(m in method_v){
    if(m=='mds'){
      z_embed <- read.csv('result/Evaluation/sim_rev_1-N200-fmds-0.00-Z.csv')
      
      res <- rp_eval2(dm=as.matrix(sim_data_dist), y_orig=sim_y_mat, z_emb=z_embed)
      c <- calcContinuityFromDist(distReference=sim_data_dist, distLowDim=dist(z_embed), kTM=k)
      t <- calcTrustworthinessFromDist(distReference = sim_data_dist, distLowDim = dist(z_embed), kTM = k)
      s <- get_stress(x_dist=sim_data_dist, z_dist=dist(z_embed))
      corr_dist <- get_pearson_corr(x_dist=sim_data_dist, z_dist=dist(z_embed))
      
      eval_df[nrow(eval_df)+1,] <- list(m, p, k, c, t, s, corr_dist, res$q, res$rho)
      next
    }
    params_v <- if(m=='fmds') params1_v else params2_v
    for(p in params_v){
      z_embed <- 
        if(m=='fmds'){
          read.csv(sprintf('result/Evaluation/sim_rev_1-N200-%s-%.2f-Z.csv', m, p))
        } else read.csv(sprintf('result/Evaluation/sim_rev_1-%s-%02d-Z.csv', m, p))
      
      res <- rp_eval2(dm=as.matrix(sim_data_dist), y_orig=sim_y_mat, z_emb=z_embed)
      c <- calcContinuityFromDist(distReference=sim_data_dist, distLowDim=dist(z_embed), kTM=k)
      t <- calcTrustworthinessFromDist(distReference = sim_data_dist, distLowDim = dist(z_embed), kTM = k)
      s <- get_stress(x_dist=sim_data_dist, z_dist=dist(z_embed))
      corr_dist <- get_pearson_corr(x_dist=sim_data_dist, z_dist=dist(z_embed))
      
      eval_df[nrow(eval_df)+1,] <- list(m, p, k, c, t, s, corr_dist, res$q, res$rho)
    }
  }
}


# Panel: Continuity v. Trustworthiness, k = 14
ggplot(data=eval_df[eval_df$k==14,]) +
  geom_point(aes(x=continuity, y=trustworthiness, shape=method, color=method, 
                 fill=as.factor(hyperparameter)), size=2.5, stroke=0.5) +
  scale_shape_manual(values=c(21,22,4,23,24,25)) + 
  # ('fmds', 'iso', 'mds', 'tsne', 'umap_s', 'umap_u')
  scale_color_manual(values=c('red','grey30','grey30','grey30','grey30','grey30')) +
  scale_fill_manual(values=c('grey99', 'grey75', 'grey25', 'grey99', 'grey75',
                             'grey25')) +
  # scale_y_continuous(breaks = seq(0.85, 1, 0.05), limits = c(0.85, 1)) +
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

ggsave('figures/fig4A_rev.pdf', width=2.41, height=1.9, units='in')


# Panel: Continuity v. Trustworthiness, k = 150
ggplot(data=eval_df[eval_df$k==150,]) +
  geom_point(aes(x=continuity, y=trustworthiness, shape=method, color=method, 
                 fill=as.factor(hyperparameter)), size=2.5, stroke=0.5) +
  scale_shape_manual(values=c(21,22,4,23,24,25)) + 
  # ('fmds', 'iso', 'mds', 'tsne', 'umap_s', 'umap_u')
  scale_color_manual(values=c('red','grey30','grey30','grey30','grey30','grey30',
                              'grey30','grey30')) +
  scale_fill_manual(values=c('grey99', 'grey75', 'grey25', 'grey99', 'grey75',
                             'grey25')) +
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

ggsave('figures/fig4B_rev.pdf', width=2.41, height=1.9, units='in')


# Panel C: stress v. corr_dist
ggplot(data=eval_df) +
  geom_point(aes(x=stress, y=corr_dist, shape=method, color=method, 
                 fill=as.factor(hyperparameter)), size=2.5, stroke=0.5) +
  scale_shape_manual(values=c(21,22,4,23,24,25)) + 
    # ('fmds', 'iso', 'mds', 'tsne', 'umap_s', 'umap_u')
  scale_color_manual(values=c('red','grey30','grey30','grey30','grey30','grey30',
                              'grey30','grey30')) +
  scale_fill_manual(values=c('grey99', 'grey75', 'grey25', 'grey99', 'grey75', 
                             'grey25')) +
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

ggsave('figures/fig4C_rev.pdf', width=2, height=1.9, units='in')


# Panel D: p_ratio v. corr_f
ggplot(data=eval_df) +
  geom_point(aes(x=corr_f, y=p_ratio, shape=method, color=method, 
                 fill=as.factor(hyperparameter)), size=2.5, stroke=0.5) +
  scale_shape_manual(values=c(21,22,4,23,24,25)) + 
    # ('fmds', 'iso', 'mds', 'tsne', 'umap_s', 'umap_u')
  scale_color_manual(values=c('red','grey30','grey30','grey30','grey30','grey30',
                              'grey30','grey30')) +
  scale_fill_manual(values=c('grey99', 'grey75', 'grey25', 'grey99', 'grey75', 
                             'grey25')) +
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

ggsave('figures/fig4D_rev.pdf', width=2.85, height=1.9, units='in')