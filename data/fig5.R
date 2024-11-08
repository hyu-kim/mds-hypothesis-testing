library(ggplot2)
library('dplyr')

dataset_v <- c('alga_2', 'sim')
iter_df <- data.frame(matrix(nrow=0, ncol=4))
colnames(iter_df) <- c('dataset', 'replicate', 'lambda', 'iteration')
for(d in dataset_v){
  for(r in 1:3){
    if(d=='alga_2' & r>1) next
    d2 <- ifelse(d=='alga_2', d, paste('sim',r, sep='_'))
    N <- 0.5 * length(list.files(sprintf('result/HyperparameterStudy/%s', d2)))
    N <- as.numeric(N)
    for(l in (4:20)/20){
      if(l==0.15 & d=='sim') {
        n_iter <- nrow(read.csv(sprintf('result/HyperparameterStudy/%s_200/%s-fmds-%.2f-log.csv', d2, d2, l)))
      }
      else n_iter <- nrow(read.csv(sprintf('result/HyperparameterStudy/%s/%s-fmds-%.2f-log.csv', d2, d2, l)))
      iter_df <- rbind(iter_df, data.frame(dataset=d, replicate=r, lambda=l, iteration=n_iter))
    }
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
  scale_x_continuous(breaks=seq(5)/5) +
  scale_y_continuous(trans='log10') +
  scale_shape_manual(values=c(1,2), labels=c("Algal", "Simulated")) + 
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
        legend.position = "right",
        axis.title.x = element_text(size=8, colour='black'),
        axis.title.y = element_text(size=8, colour='black'),
        axis.text.x = element_text(size=8, colour='black'),
        axis.text.y = element_text(size=8, colour='black'),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth=0.25, colour = 'black')
  )

ggsave('figures/fig5.pdf', width=4, height=2, units='in')