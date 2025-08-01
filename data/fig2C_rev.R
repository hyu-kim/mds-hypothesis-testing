library(ggplot2)
library('dplyr')

dataset_v <- c('sim')
iter_df <- data.frame(matrix(nrow=0, ncol=4))
colnames(iter_df) <- c('dataset', 'replicate', 'lambda', 'iteration')

for(r in 1:3){
  for(l in (2:20)/20){
    n_iter <- nrow(read.csv(sprintf('result/HyperparameterStudy/sim_rev/sim_rev_%d-fmds-%.2f-log.csv', r, l)))
    iter_df <- rbind(iter_df, data.frame(dataset='sim_rev', replicate=r, lambda=l, iteration=n_iter))
  }
}


iter_df_stat <- iter_df %>%
  group_by(dataset, lambda) %>%
  summarize(iter_mean = mean(iteration), iter_std = sd(iteration), n = n())


ggplot(data=iter_df_stat) +
  geom_point(aes(x=lambda, y=iter_mean-1, shape=dataset)) +
  geom_errorbar(aes(x=lambda, ymin=iter_mean-1-iter_std, ymax=iter_mean-1+iter_std),
                width=0.03, color='black', size=0.25) +
  # geom_line(aes(x=lambda, y=iteration, color=dataset), size=0.25, linetype='solid') +
  scale_x_continuous(limits = c(0,1), breaks=(0:5)/5) +
  scale_y_continuous(trans='log10') +
  scale_shape_manual(values=1, labels="Simulated") +
  labs(x ="Hyperparameter", y = "Iteration") +
  # scale_color_manual(values=c('black','black')) +
  theme(strip.background = element_rect(fill=NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(fill = "transparent", color = 'black', size=0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=8),
        legend.position = "none",
        axis.title.x = element_text(size=8, colour='black'),
        axis.title.y = element_text(size=8, colour='black'),
        axis.text.x = element_text(size=8, colour='black'),
        axis.text.y = element_text(size=8, colour='black'),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth=0.25, colour = 'black')
  )

ggsave('figures/fig2C_rev.pdf', width=4, height=2, units='in')