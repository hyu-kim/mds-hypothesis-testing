library(ggplot2)
library(dplyr)
library(ggnewscale)
library(tidyr)
library(vegan)

x_mat <- readRDS('result/Multiclass/sim4d_rev-data.rds')

# transform to log scale for visibility
log_mat <- x_mat
log_mat[x_mat == 0] <- NA
log_mat[!is.na(log_mat)] <- log(log_mat[!is.na(log_mat)])

# change to long type for ggplot
mat_long <- as.data.frame(log_mat) %>%
  mutate(row = rownames(log_mat)) %>%
  pivot_longer(-row, names_to = "col", values_to = "value")
last_two <- substr(mat_long$col, nchar(mat_long$col) - 1, nchar(mat_long$col))
mat_long$"group" <- ifelse(last_two == '-2', 2, ifelse(last_two == '-3', 3, 1))
mat_long$col <- factor(mat_long$col, levels = colnames(log_mat))


## Panel A
ggplot() +
  # Heatmap for group A
  geom_tile(data = filter(mat_long, group == "1"),
            aes(row, col, fill = value)) +
  scale_fill_gradient(low = "white", high = "red", name = "Group 1", na.value = "white", limits = c(-20, 0),
                      labels = c(expression(10^{-20}), expression(10^{-15}), expression(10^{-10}), expression(10^{-5}), expression(10^{0})),
                      guide = guide_colorbar(frame.colour = "black", frame.linewidth = 0.1, ticks.colour = "black", ticks.linewidth = 0.1)) +
  # Start new fill scale for group B
  new_scale_fill() +
  # Heatmap for group B
  geom_tile(data = filter(mat_long, group == "2"),
            aes(row, col, fill = value)) +
  scale_fill_gradient(low = "white", high = "blue", name = "Group 2", na.value = "white", limits = c(-20, 0),
                      labels = c(expression(10^{-20}), expression(10^{-15}), expression(10^{-10}), expression(10^{-5}), expression(10^{0})),
                      guide = guide_colorbar(frame.colour = "black", frame.linewidth = 0.1, ticks.colour = "black", ticks.linewidth = 0.1)) +
  new_scale_fill() +
  # Heatmap for group C
  geom_tile(data = filter(mat_long, group == "3"),
            aes(row, col, fill = value)) +
  scale_fill_gradient(low = "white", high = "green4", name = "Group 3", na.value = "white", limits = c(-20, 0),
                      labels = c(expression(10^{-20}), expression(10^{-15}), expression(10^{-10}), expression(10^{-5}), expression(10^{0})),
                      guide = guide_colorbar(frame.colour = "black", frame.linewidth = 0.1, ticks.colour = "black", ticks.linewidth = 0.1)) +
  labs(x = "Microbial taxon (species)", y = "Sample", fill = "Group") +
  theme(strip.background = element_rect(fill=NA),
        panel.background = element_rect(color = 'black', linewidth = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        # panel.border = element_rect(fill = "transparent", color = 'black', size=0.5),
        legend.position = "bottom",
        legend.key.size = unit(0.25, "in"),
        legend.title = element_text(size = 10, hjust = 0.5, vjust = 0.9),
        legend.text = element_text(size = 10),
        # legend.background = element_rect(fill="grey"),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_line(color = "grey25", linewidth = 0.5),
        axis.ticks = element_line(color = "grey25", linewidth = 0.1))

ggsave('figures/Fig_S8A_rev.pdf', width=7, height=5, units='in')


## Panel B
mat_long_1 <- as.numeric(filter(mat_long, group == "1")$value)
mat_long_1 <- mat_long_1[!is.na(mat_long_1)]

mat_long_2 <- as.numeric(filter(mat_long, group == "2")$value)
mat_long_2 <- mat_long_2[!is.na(mat_long_2)]

mat_long_3 <- as.numeric(filter(mat_long, group == "3")$value)
mat_long_3 <- mat_long_3[!is.na(mat_long_3)]

ggplot(mat_long[!is.na(mat_long$value),], aes(x = value, fill = as.factor(group))) +
  geom_histogram(aes(y = after_stat(density)), position = "identity", alpha = 0.5, bins = 30) +
  scale_fill_manual(values = c("red", "blue", "green4")) +
  labs(x = "Log relative abundance", y = "Density", fill = "Group") +
  theme(strip.background = element_rect(fill=NA),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        # panel.border = element_rect(fill = "transparent", color = 'black', size=0.5),
        legend.position = "bottom",
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        # legend.background = element_rect(fill="grey"),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.5))

ggsave('figures/Fig_S8B_rev.pdf', width=4, height=3, units='in')


## Panel C
dist <- vegdist(t(x_mat), method="bray")
pcoa <- cmdscale(dist, eig = TRUE)
pcoa_df <- as.data.frame(pcoa$points)
group <- factor(rep(c(1, 2, 3), each = 25))
pcoa_df$group <- group
dist2 <- dist(pcoa$points[,1:2], method='euclidean')
result <- adonis2(dist ~ group)
px <- result$`Pr(>F)`[1]
result2 <- adonis2(dist2 ~ group)
pz <- result2$`Pr(>F)`[1]

ggplot(pcoa_df) +
  geom_point(aes(x = V1, y = V2, color = group), size = 2, alpha = 0.7) +
  scale_color_manual(values = c("1" = "red", "2" = "blue", "3" = "green4")) +
  labs(x = "V1", y = "V2", color = "Group") +
  theme(strip.background = element_rect(fill=NA),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        # panel.border = element_rect(fill = "transparent", color = 'black', size=0.5),
        legend.position = "bottom",
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        # legend.background = element_rect(fill="grey"),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.5))

ggsave('figures/Fig_S8C_rev.pdf', width=2.7, height=3, units='in')