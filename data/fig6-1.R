library(ggplot2)

method_v <- c('fmds', 'mds')
params_v <- c(2,4)/4

vis_df <- data.frame(matrix(ncol=5, nrow=0))
colnames(vis_df) <- c('method', 'hyperparameter', 'y', 'X1', 'X2')

sim_data <- as.matrix(read.csv('result/Evaluation/sim_1-data.csv'))
y_sim <- as.matrix(read.csv('result/Evaluation/sim_1-Y.csv'))

for(m in method_v){
  if(m=='mds'){
    z_embed <- read.csv('result/Evaluation/sim_1-mds-Z.csv')
    vis_df <- rbind(vis_df, data.frame(method = m, hyperparameter = 0, 
                                       y = as.factor(y_sim), 
                                       X1 = z_embed[,1], X2 = z_embed[,2]))
    next
  }
  for(p in params_v){
    z_embed <- read.csv(sprintf('result/Evaluation/sim_1-%s-%.2f-Z.csv', m, p))
    vis_df <- rbind(vis_df, data.frame(method = m, hyperparameter = p, 
                                       y = as.factor(y_sim), 
                                       X1 = z_embed[,1], X2 = z_embed[,2]))
  }
}


# plot
ggplot(data=vis_df, aes(x=X1, y=X2, shape=y, linetype=y)) +
  geom_point(aes(fill=y), size=1.75, stroke=0, alpha=0.8) +
  stat_ellipse(aes(fill=y), geom='polygon', lwd=0.2, level=0.68, alpha=0.25) +
  facet_wrap(facets = vars(hyperparameter), ncol=4) +
  scale_shape_manual(values=c(21,24)) +
  scale_fill_manual(values=c('red', 'blue')) +
  scale_linetype_manual(values=c('longdash','twodash')) +
  # scale_x_continuous(limits = c(-5, 5), breaks = seq(-4, 4, 2)) +
  # scale_y_continuous(limits = c(-4, 4), breaks = seq(-4, 4, 2)) +
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
        axis.ticks = element_line(linewidth=0.25, colour = 'black'),
        strip.text.x = element_blank()
        )

ggsave('figures/fig6A.pdf', width=5.1, height=1.78, units='in')