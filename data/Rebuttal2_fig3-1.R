library(ggplot2)
library(tidyr)
library(lemon)
library(cowplot)
library(patchwork)
library(ggh4x)
library(dplyr)  # Added for data manipulation

## A. without stopping rule
log_all_df <- data.frame(matrix(nrow=0, ncol=6))
colnames(log_all_df) <- c('dataset', 'size', 'lambda', 'epoch', 'p_z', 'p_0')

for(r in c(1)){
  for(N in c(100)){
    for(l in c(0.2, 0.5, 0.8)){
      log_df <- read.csv(sprintf('result/Revision2_Dec2025/sim_rev_%g-N%g-fmds-%.2f-log.csv', r, N, l))
      log_all_df <- rbind(log_all_df, 
                          data.frame(dataset = sprintf('sim_rev_%g', r), 
                                     size = N, lambda = sprintf('lambda = %g', l), 
                                     epoch = log_df$epoch, 
                                     p_z = log_df$p_z, p_0 = log_df$p_0
                          ))
    }
  }
}

## Identify first convergence point for each lambda
first_convergence <- log_all_df %>%
  group_by(lambda) %>%
  mutate(within_threshold = abs(p_z - p_0) < 0.05) %>%
  filter(within_threshold) %>%
  slice(1) %>%  # Take the first occurrence
  ungroup()

ggplot(log_all_df) + 
  geom_point(aes(x=epoch, y=p_z), size = 0.5) +
  geom_line(aes(x=epoch, y=p_z), linewidth = 0.15) +
  geom_hline(aes(yintercept = p_0 - 0.05), linetype = "dashed", color = "red", linewidth = 0.3) +
  geom_hline(aes(yintercept = p_0 + 0.05), linetype = "dashed", color = "red", linewidth = 0.3) +
  # Highlight first convergence point
  geom_point(data = first_convergence, 
             aes(x = epoch, y = p_z), 
             color = "red", size = 2, shape = 21, fill = "magenta", stroke = 0.5) +
  facet_wrap2(~lambda, scales = "free", nrow = 1) +
  labs(x = "Epoch", y = bquote(p[z])) +
  scale_y_continuous(limits = c(0, 1), breaks=seq(0,1,0.25)) +
  theme(strip.background = element_rect(fill=NA),
        # strip.text = element_blank(),
        text = element_text(size = 8),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # plot.background = element_blank(),
        panel.border = element_rect(fill = "transparent", color = 'black', size=0.5),
        legend.position = "Bottom",
        axis.text.x = element_text(size=8, colour='black'),
        axis.text.y = element_text(size=8, colour='black'),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth=0.25, colour = 'black')
  )

ggsave('figures/Rebuttal2_Fig3-2a.png', width=6, height=2.5, units='in')


## B. with stopping rule
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
                                     size = sprintf("N = %g", N), lambda = l, 
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
  facet_grid(dataset~size) +
  scale_x_continuous(expand = c(0.003, 0)) +
  scale_y_continuous(expand = c(0.005, 0)) +
  labs(x = "Epoch", y = bquote(lambda)) +
  scale_fill_gradient(low = "white", high = "black", limits = c(0,1), 
                      breaks = c(0, 0.5, 1), name = bquote(p[z]),
                      guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", linewidth = 0.25)) +
  theme(text = element_text(size = 8),
        panel.background = element_rect(fill = 'white'),
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
        legend.position="right",
        legend.text=element_text(size=8),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.15, "in")
  )

ggsave('figures/Rebuttal2_Fig3-2b.png', width=6, height=4, units='in')