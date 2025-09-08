library(ggplot2)
library(tidyr)
library(lemon)

v_lambda = c(0, 1, 2, 4, 10)/10
v_linestyle = c('dotted', 'dashed', 'solid', 'blank', 'blank')
v_linesize = c(0.5, 0.4, 0.2, 0.1, 0.1)
# v_color = c('#f8766d','#7cae00','#00bfc4','black','black')

log_all_df <- data.frame(matrix(nrow=0, ncol=7))
colnames(log_all_df) <- c('dataset', 'size', 'lambda', 'epoch', 'p_z', 'p_0', 'tile_size')

for(r in c(1:3)){
  for(N in c(50, 100, 200, 500)){
    for(l in c((0:4)/50, (2:20)/20)){
        # if(r==3 & N==50 & l<0.05){next}
      log_df <- 
        read.csv(sprintf('result/ScalingStudy/sim_rev_%g/sim_rev_%g-N%d-fmds-%.2f-log.csv',
                         r, r, N, l))
      log_all_df <- rbind(log_all_df, 
                          data.frame(dataset = sprintf('replicate %g', r), 
                                     size = N, lambda = l, 
                                     epoch = log_df$epoch, p_z = log_df$p_z, 
                                     p_0 = log_df$p_0, 
                                     tile_size = ifelse(l<0.1, 0.02, 0.05)))
    }
  }
}

# adjust tiles at l = 0.1
log_all_df$tile_size[log_all_df$lambda==0.1] <- 0.035
log_all_df$lambda[log_all_df$lambda==0.1] <- 0.107

ggplot(data = log_all_df) +
  geom_tile(aes(x = epoch, y = lambda, fill = p_z, height = tile_size)) +
  facet_grid(size~dataset) +
  scale_x_continuous(expand = c(0.003, 0)) +
  scale_y_continuous(expand = c(0.005, 0)) +
  scale_fill_gradient(low = "white", high = "black", limits = c(0,1), 
                      breaks = c(0, 0.5, 1),
                      guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", linewidth = 0.25)) +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        axis.line = element_line(colour = "black", size = 0.25),
        axis.title.x = element_text(size=8, colour='black'),
        axis.title.y = element_text(size=8, colour='black'),
        axis.text.x = element_text(size=8, colour='black'),
        axis.text.y = element_text(size=8, colour='black'),
        axis.text.x.top = element_blank(),
        axis.text.y.right = element_blank(),
        axis.title.x.top = element_blank(),
        axis.title.y.right = element_blank(),
        legend.position="bottom",
        legend.text=element_text(size=8),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.15, "in")
        )

ggsave('figures/Fig_S2_rev.pdf', width=6, height=7, units='in')