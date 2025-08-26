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

ggsave('figures/Fig_S1_rev.pdf', width=6, height=7, units='in')


# ## Panel A
# # import
# v_lambda_1 = (0:3)/20 # lambdas for 200-iterated data
# v_lambda_2 = c(0,1,2,5,10)/10 # lambdas for 50-iterated data
# ratio_v = c(100, 1, 5)
# lambda = 0 # 0, 0.2, 0.4
# display_ratio = 0.2 # 0.2, 0.2
# df_log <- read.csv(paste('result/HyperparameterStudy/sim_1',
#                          # ifelse(lambda>=0.10, '50', '200'),
#                          '/sim_1-fmds-', sprintf('%.2f',lambda), '-log.csv', sep=''))
# 
# ggplot(df_log) + 
#   # geom_point(aes(x=epoch, y=obj_mds), size=0.75, shape=1) + # iter_50 only
#   # geom_point(aes(x=epoch, y=obj_mds), inherit.aes = FALSE, size=1, shape=4) +
#   # geom_point(aes(x=epoch, y=obj_confr*0.2), inherit.aes = FALSE, size=1, shape=1) +
#   geom_line(aes(x=epoch, y=obj_mds), size=0.25, linetype='dashed') +
#   geom_line(aes(x=epoch, y=obj_confr*0.2), inherit.aes = FALSE, size=0.25, linetype='solid') + # iter_50 only
#   scale_y_continuous(sec.axis = sec_axis(~./0.2, name="Confirmatory", 
#     breaks=c(0,50,100)), limits=c(0,24), breaks = c(0,10,20)) +
#   theme(strip.background = element_rect(fill=NA),
#         panel.background = element_rect(fill = "transparent", color = NA),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         plot.background = element_blank(),
#         panel.border = element_rect(fill = "transparent", color = 'black', size=0.5),
#         legend.position = "None",
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text.x = element_text(size=8, colour='black'),
#         axis.text.y = element_text(size=8, colour='black'),
#         axis.line = element_blank(),
#         axis.ticks = element_line(linewidth=0.25, colour = 'black')
#   )
# 
# ggsave(paste('figures/fig1A_', sprintf('%.2f', lambda), '_obj.pdf', sep=''), 
#        width=2.6, height=1.4, units='in')