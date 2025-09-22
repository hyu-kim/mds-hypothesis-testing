library(vegan)
library(ggplot2)
library(dplyr)
source('fig_util.R')
source('../source/dreval/continuity.R')
source('../source/dreval/trustworthiness.R')

method_v <- c('fmds', 'mds', 'umap_s', 'umap_u', 'tsne', 'iso')
params1_v <- c(0.2, 0.4, 0.8)
params2_v <- c(5, 10, 50)

eval_df <- data.frame(matrix(ncol=11, nrow=0))
colnames(eval_df) <- c('N', 'replicate', 'method', 'hyperparameter', 'k', 'continuity', 
                       'trustworthiness', 'stress', 'corr_dist', 'p_ratio', 'corr_f')

sim_y_mat <- as.matrix(read.csv(sprintf('result/ScalingStudy/sim_rev-N200-Y.csv')))

for(r in seq(3)){
  sim_data <- t(readRDS(sprintf('result/ScalingStudy/sim_rev_%g/sim_rev_%g-N200-data.Rds', r, r)))
  rownames(sim_data) <- c(1:200)
  sim_data_dist <- vegdist(sim_data, method="bray")
  
  for(k in c(14,150)){ # local, global
  message('Replicate: ', r, ', k: ', k)
    for(m in method_v){
      if(m=='mds'){
        z_embed <- read.csv(sprintf('result/Evaluation/SimRev%g/sim_rev_%g-N200-fmds-0.00-Z.csv',r, r))
        
        res <- rp_eval2(dm=as.matrix(sim_data_dist), y_orig=sim_y_mat, z_emb=z_embed)
        c <- calcContinuityFromDist(distReference=sim_data_dist, distLowDim=dist(z_embed), kTM=k)
        t <- calcTrustworthinessFromDist(distReference = sim_data_dist, distLowDim = dist(z_embed), kTM = k)
        s <- get_stress(x_dist=sim_data_dist, z_dist=dist(z_embed))
        corr_dist <- get_pearson_corr(x_dist=sim_data_dist, z_dist=dist(z_embed))
        
        eval_df[nrow(eval_df)+1,] <- list(200, r, m, p, k, c, t, s, corr_dist, res$q, res$rho)
        next
      }
      params_v <- if(m=='fmds') params1_v else params2_v
      for(p in params_v){
        z_embed <- 
          if(m=='fmds'){
            read.csv(sprintf('result/Evaluation/SimRev%g/sim_rev_%g-N200-%s-%.2f-Z.csv', r, r, m, p))
          } else read.csv(sprintf('result/Evaluation/SimRev%g/sim_rev_%g-%s-%02d-Z.csv', r, r, m, p))
        
        res <- rp_eval2(dm=as.matrix(sim_data_dist), y_orig=sim_y_mat, z_emb=z_embed)
        c <- calcContinuityFromDist(distReference=sim_data_dist, distLowDim=dist(z_embed), kTM=k)
        t <- calcTrustworthinessFromDist(distReference = sim_data_dist, distLowDim = dist(z_embed), kTM = k)
        s <- get_stress(x_dist=sim_data_dist, z_dist=dist(z_embed))
        corr_dist <- get_pearson_corr(x_dist=sim_data_dist, z_dist=dist(z_embed))
        
        eval_df[nrow(eval_df)+1,] <- list(200, r, m, p, k, c, t, s, corr_dist, res$q, res$rho)
      }
    }
  }
}

# save or load Rds
saveRDS(eval_df, "result/fig4_rev_evaldf2.Rds")
eval_df <- readRDS("result/fig4_rev_evaldf.Rds")

eval_df_stat <- eval_df %>%
  group_by(method, hyperparameter, k) %>%
  summarize(cont_mu = mean(continuity), cont_sd = sd(continuity),
            trus_mu = mean(trustworthiness), trus_sd = sd(trustworthiness),
            stre_mu = mean(stress), stre_sd = sd(stress),
            corr_d_mu = mean(corr_dist), corr_d_sd = sd(corr_dist),
            p_rat_mu = mean(p_ratio), p_rat_sd = sd(p_ratio),
            corr_f_mu = mean(corr_f), corr_f_sd = sd(corr_f),
            n = n())

eval_df_viz <- eval_df_stat[eval_df_stat$hyperparameter!=0.4 & eval_df_stat$hyperparameter!=10,]


# Panel A: Continuity v. Trustworthiness, k = 14
ggplot(data = eval_df_viz[eval_df_viz$k==14,]) +
  geom_point(aes(x=cont_mu, y=trus_mu, shape=method, color=method, 
                 fill=as.factor(hyperparameter)), size=2.5, stroke=0.5) +
  geom_errorbar(aes(x=cont_mu, ymin=trus_mu-trus_sd, ymax=trus_mu+trus_sd),
                width=0.005, color='grey25', size=0.25) +
  geom_errorbarh(aes(xmin=cont_mu-cont_sd, xmax=cont_mu+cont_sd, y=trus_mu),
                 height=0.005, color='grey25', size=0.25) +
  scale_shape_manual(values=c(21,22,4,23,24,25)) + 
  # ('fmds', 'iso', 'mds', 'tsne', 'umap_s', 'umap_u')
  scale_color_manual(values=c('red','grey20','grey20','grey20','grey20','grey20')) +
  scale_fill_manual(values=c('grey99', 'grey50', 'grey99', 'grey40')) +
  scale_x_continuous(breaks=c(0.6,0.7,0.8)) +
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
ggplot(data = eval_df_viz[eval_df_viz$k==150,]) +
  geom_point(aes(x=cont_mu, y=trus_mu, shape=method, color=method, 
                 fill=as.factor(hyperparameter)), size=2.5, stroke=0.5) +
  geom_errorbar(aes(x=cont_mu, ymin=trus_mu-trus_sd, ymax=trus_mu+trus_sd),
                width=0.005, color='grey25', size=0.25) +
  geom_errorbarh(aes(xmin=cont_mu-cont_sd, xmax=cont_mu+cont_sd, y=trus_mu),
                 height=0.005, color='grey25', size=0.25) +
  scale_shape_manual(values=c(21,22,4,23,24,25)) + 
  # ('fmds', 'iso', 'mds', 'tsne', 'umap_s', 'umap_u')
  scale_color_manual(values=c('red','grey20','grey20','grey20','grey20','grey20')) +
  scale_fill_manual(values=c('transparent', 'grey50', 'transparent', 'grey40')) +
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
ggplot(data = eval_df_viz[eval_df_viz$k==150,]) +
  geom_point(aes(x=stre_mu, y=corr_d_mu, shape=method, color=method, 
                 fill=as.factor(hyperparameter)), size=2.5, stroke=0.5) +
  geom_errorbar(aes(x=stre_mu, ymin=corr_d_mu-corr_d_sd, ymax=corr_d_mu+corr_d_sd),
                width=0.015, color='grey25', size=0.25) +
  geom_errorbarh(aes(xmin=stre_mu-stre_sd, xmax=stre_mu+stre_sd, y=corr_d_mu),
                height=0.015, color='grey25', size=0.25) +
  scale_shape_manual(values=c(21,22,4,23,24,25)) + 
    # ('fmds', 'iso', 'mds', 'tsne', 'umap_s', 'umap_u')
  scale_color_manual(values=c('red','grey20','grey20','grey20','grey20','grey20')) +
  scale_fill_manual(values=c('transparent', 'grey50', 'transparent', 'grey40')) +
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
ggplot(data = eval_df_viz[eval_df_viz$k==14,]) + # k doesn't matter
  geom_point(aes(x=corr_f_mu, y=p_rat_mu, shape=method, color=method, 
                 fill=as.factor(hyperparameter)), size=2.5, stroke=0.5) +
  geom_errorbar(aes(x=corr_f_mu, ymin=p_rat_mu-p_rat_sd, ymax=p_rat_mu+p_rat_sd),
                width=0.01, color='grey25', size=0.25) +
  geom_errorbarh(aes(xmin=corr_f_mu-corr_f_sd, xmax=corr_f_mu+corr_f_sd, y=p_rat_mu),
                height=0.025, color='grey25', size=0.25) +
  scale_shape_manual(values=c(21,22,4,23,24,25)) + 
    # ('fmds', 'iso', 'mds', 'tsne', 'umap_s', 'umap_u')
  scale_color_manual(values=c('red','grey20','grey20','grey20','grey20','grey20')) +
  scale_fill_manual(values=c('transparent', 'grey50', 'transparent', 'grey40')) +
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