library(ggplot2)
library(ggh4x)

## Panel A
N <- 200 # exemplary data where N ~ d
log_all_df <- data.frame(matrix(nrow=0, ncol=8))
colnames(log_all_df) <- c('dataset', 'size', 'lambda', 'epoch', 'obj_mds', 'obj_confr', 'p_z', 'p_0')

for(l in c(0, 0.04, 0.2, 0.8)){
  log_df <- read.csv(sprintf('result/ScalingStudy/sim_rev_1/sim_rev_1-N200-fmds-%.2f-log.csv', l))
  log_all_df <- rbind(log_all_df, 
                      data.frame(dataset = 'sim_rev_1', 
                                 size = N, lambda = l, 
                                 epoch = log_df$epoch, 
                                 obj_mds = log_df$obj_mds, obj_confr = log_df$obj_confr, 
                                 p_z = log_df$p_z, p_0 = log_df$p_0
                                 ))
}

# filter rows for better display
log_all_df <- log_all_df[log_all_df$epoch<31,]
log_all_df <- log_all_df[log_all_df$lambda>0.1 | log_all_df$epoch%%2!=1,]

ggplot(log_all_df) + 
  geom_point(aes(x=epoch, y=obj_mds), shape=1, size=1.5) +
  geom_point(aes(x=epoch, y=obj_confr), shape=2, size=1.5) +
  facet_wrap(~lambda, nrow = 1, scales = "free") +
  scale_y_continuous(limits = c(0, 600), breaks=seq(0,600,200)) +
  theme(strip.background = element_rect(fill=NA),
        strip.text = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(fill = "transparent", color = 'black', size=0.5),
        legend.position = "Bottom",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=8, colour='black'),
        # axis.text.y = element_text(size=8, colour='black'),
        axis.text.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth=0.25, colour = 'black')
  )

ggsave('figures/fig2A_rev.pdf', width=4.8, height=1.35, units='in')



## Panel B
log_all_df <- data.frame(matrix(nrow=0, ncol=6))
colnames(log_all_df) <- c('dataset', 'size', 'lambda', 'epoch', 'p_z', 'p_0')

for(r in c(1)){
  for(N in c(50,100,200,500)){
    for(l in c(0.04, 0.2, 0.8)){
      log_df <- read.csv(sprintf('result/ScalingStudy/sim_rev_%g/sim_rev_%g-N%g-fmds-%.2f-log.csv', r, r, N, l))
      log_all_df <- rbind(log_all_df, 
                          data.frame(dataset = sprintf('sim_rev_%g', r), 
                                     size = N, lambda = l, 
                                     epoch = log_df$epoch, 
                                     p_z = log_df$p_z, p_0 = log_df$p_0
                          ))
    }
  }
}

# filter rows for better display
log_all_df <- log_all_df[log_all_df$epoch<32,]
log_all_df <- log_all_df[log_all_df$lambda!=0.04 | log_all_df$epoch%%2!=1,]

ggplot(log_all_df) + 
  geom_point(aes(x=epoch, y=p_z, shape=as.factor(lambda)), size=1.5) +
  facet_wrap2(~size, scales = "free", nrow = 1) +
  # facetted_pos_scales(
  #   x = list(
  #     scale_x_continuous(limits = c(0, 30)),
  #     scale_x_continuous(limits = c(0, 20)),
  #     scale_x_continuous(limits = c(0, 32)),
  #     scale_x_continuous(limits = c(0, 42))
  #   )
  # ) +
  scale_y_continuous(limits = c(0, 1), breaks=seq(0,1,0.25)) +
  scale_shape_manual(values=c(21,22,24)) +
  theme(strip.background = element_rect(fill=NA),
        strip.text = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(fill = "transparent", color = 'black', size=0.5),
        legend.position = "Bottom",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=8, colour='black'),
        # axis.text.y = element_text(size=8, colour='black'),
        axis.text.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth=0.25, colour = 'black')
  )

ggsave('figures/fig2B_rev.pdf', width=4.8, height=1.35, units='in')



## Panel C
iter_df <- data.frame(matrix(nrow=0, ncol=4))
colnames(iter_df) <- c('dataset', 'size', 'lambda', 'iteration')

for(r in c(1:3)){
  for(N in c(50, 100, 200, 500)){
    for(l in c((1:10)/10)){
      if(r==3 & N==50 & l<0.05){next} # did not converge
      n_iter <- nrow(read.csv(sprintf('result/ScalingStudy/sim_rev_%g/sim_rev_%g-N%d-fmds-%.2f-log.csv',
                                      r, r, N, l)))
      iter_df <- rbind(iter_df, 
                       data.frame(dataset=sprintf('replicate %g', r), 
                                  size=N, lambda=l, iteration=n_iter))
    }
  }
}

iter_df <- iter_df[iter_df$iteration!=51,]

iter_df_stat <- iter_df %>%
  group_by(lambda, size) %>%
  summarize(iter_mean = mean(iteration), iter_std = sd(iteration), n = n())

iter_df_stat <- iter_df_stat[iter_df_stat$n==3,]

ggplot(data=iter_df_stat) +
  geom_point(aes(x=lambda, y=iter_mean-1), shape=1, size=1.5) +
  geom_errorbar(aes(x=lambda, ymin=iter_mean-1-iter_std, ymax=iter_mean-1+iter_std),
                width=0.05, color='black', size=0.25) +
  facet_wrap2(~size, axes="all", nrow=1) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_x_continuous(limits = c(0.05,1.1), breaks=(0:2)/2) +
  # scale_y_continuous(limits = c(1, 45), labels = scales::label_math(10^.x)) +
  annotation_logticks(sides = "l", outside = TRUE, linewidth = 0.25,
                      short = unit(0.7, "mm"),  # longer minor ticks
                      mid = unit(1, "mm"),      # medium ticks
                      long = unit(1.5, "mm")) +      # longest major ticks)
  coord_cartesian(clip = "off") +
  labs(x ="Hyperparameter", y = "No. iteration") +
  # scale_color_manual(values=c('black','black')) +
  theme(strip.background = element_rect(fill=NA),
        strip.text = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(fill = "transparent", color = 'black', size=0.25),
        legend.title=element_blank(), 
        legend.text=element_text(size=8),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=8, colour='black'),
        axis.text.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks.x = element_line(linewidth=0.25, colour = 'black'),
        axis.ticks.y = element_blank()
  )

ggsave('figures/fig2C_rev.pdf', width=4.8, height=1.35, units='in')