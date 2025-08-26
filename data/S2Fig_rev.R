library(ggplot2)
library('dplyr')

iter_df <- data.frame(matrix(nrow=0, ncol=4))
colnames(iter_df) <- c('dataset', 'size', 'lambda', 'iteration')

for(r in c(1:3)){
  for(N in c(50, 100, 200, 500)){
    for(l in c((0:4)/50, (2:20)/20)){
      if(r==3 & N==50 & l<0.05){next} # did not converge
      n_iter <- nrow(read.csv(sprintf('result/ScalingStudy/sim_rev_%g/sim_rev_%g-N%d-fmds-%.2f-log.csv',
                                      r, r, N, l)))
      iter_df <- rbind(iter_df, 
                       data.frame(dataset=sprintf('replicate %g', r), 
                                  size=N, lambda=l, iteration=n_iter))
    }
  }
}

# filter out if n_iter = 51 (max)
print(iter_df[iter_df$iteration==51,])
iter_df <- iter_df[iter_df$iteration!=51,]


ggplot(data=iter_df) +
  geom_point(aes(x=lambda, y=iteration-1, fill=as.factor(size)), shape = 21, color = 'black', size=1.5) +
  # geom_errorbar(aes(x=lambda, ymin=iter_mean-1-iter_std, ymax=iter_mean-1+iter_std),
  #               width=0.03, color='black', size=0.25) +
  geom_line(aes(x=lambda, y=iteration-1, color = as.factor(size))) +
  facet_grid(~dataset) +
  scale_x_continuous(limits = c(0,1), breaks=(0:5)/5) +
  scale_y_continuous(trans='log10') +
  # scale_shape_manual(values= c(21)) +
  # scale_fill_brewer(palette = "pal10") +
  scale_fill_manual(values=c('gold','olivedrab', 'blue', 'grey25')) +
  scale_color_manual(values=c('gold','olivedrab', 'blue', 'grey25')) +
  labs(x ="Hyperparameter", y = "Number of iterations") +
  theme(strip.background = element_rect(fill=NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(fill = "transparent", color = 'black', size=0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=8),
        legend.position = "bottom",
        axis.title.x = element_text(size=8, colour='black'),
        axis.title.y = element_text(size=8, colour='black'),
        axis.text.x = element_text(size=8, colour='black'),
        axis.text.y = element_text(size=8, colour='black'),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth=0.25, colour = 'black')
  )

ggsave('figures/Fig_S2_rev.pdf', width=6.5, height=4, units='in')