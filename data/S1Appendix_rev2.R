# loads sim_1_rev data
# permute labels, compute F ratios in Z and X
# compare F ratios ordered vs. as-is, two panels
# compare F ratios permuted once for all vs. twice for each
# highlight point where un-permuted

library(vegan)
library(scales)
source("permanova_with_config.R")
library(ggplot2)
library(cowplot)
library(patchwork)
library(ggh4x)


get_y_rand <- function(y){
  N <- 200
  a <- 2
  tbl <- table(y)
  y_rand <- rep(1,N)
  pool_v <- 1:N
  for(cl in 2:a){
    ind_rand <- sample(pool_v, tbl[cl], replace=F)
    y_rand[ind_rand] = cl
    pool_v <- setdiff(pool_v, ind_rand)
  }
  
  return(y_rand)
}

get_f_pair <- function(x_mat, y){
  d_x <- vegdist(t(x_mat), method="bray")
  pcoa <- cmdscale(d_x, eig = TRUE)
  pcoa_df <- as.data.frame(pcoa$points)
  d_z <- dist(pcoa_df, method='euclidean')
  f_x <- pseudo_F(d = d_x, trt = y)
  f_z <- pseudo_F(d = d_z, trt = y)
  
  return(list(f_x = f_x, f_z = f_z))
}

permute_pair_f <- function(n_iter, x_mat, y){
  f_paired_mat = matrix(0, nrow=n_iter+1, ncol=2) # 1st, x; 2nd, z
  f_paired <- get_f_pair(x_mat, y)
  f_paired_mat[1,1]<- f_paired$f_x
  f_paired_mat[2,2]<- f_paired$f_z
  
  for(i in 2:(n_iter+1)){
    y_rand <- get_y_rand(y)
    f_paired <- get_f_pair(x_mat, y_rand)
    f_paired_mat[i,1]<- f_paired$f_x
    f_paired_mat[i,2]<- f_paired$f_z
  }
  
  return(f_paired_mat)
}

permute_sort_f <- function(n_iter, x_mat, y){ # once for both x and z
  f_paired_mat <- permute_pair_f(n_iter, x_mat, y)  
  f_sorted_x = sort(f_paired_mat[2:(n_iter+1),1], decreasing = TRUE)
  f_sorted_z = sort(f_paired_mat[2:(n_iter+1),2], decreasing = TRUE)
  f_sorted_mat <- cbind(f_sorted_x, f_sorted_z)
  
  return(list(sorted = f_sorted_mat, paired = f_paired_mat))
}


## Panel A-C
x_mat <- readRDS('result/ScalingStudy/sim_rev_1/sim_rev_1-N200-data.Rds')
y <- read.csv('result/ScalingStudy/sim_rev-N200-Y.csv')$x

n_iter <- 999

par(mfrow = c(1, 3))

res <- permute_sort_f(n_iter, x_mat, y)
f_sorted_mat <- res$sorted
f_paired_mat <- res$paired

p1 <- ggplot() +
  geom_point(aes(x=f_sorted_mat[,1], y=f_sorted_mat[,2]), size=0.5, alpha = 0.4) +
  geom_point(aes(x=f_paired_mat[1,1], y=f_paired_mat[1,2]), col = "red", pch = 4, size = 3) +
  labs(x = bquote("Sorted F"[x]^pi), y = bquote("Sorted F"[z]^pi)) +
  theme(strip.background = element_rect(fill=NA),
        text = element_text(size = 8),
        strip.text = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(fill = "transparent", color = 'black', size=0.5),
        legend.position = "Bottom",
        axis.title.x = element_text(size=8, colour='black'),
        axis.title.y = element_text(size=8, colour='black'),
        axis.text.x = element_text(size=8, colour='black'),
        axis.text.y = element_text(size=8, colour='black'),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth=0.25, colour = 'black')
  )

p2 <- ggplot() +
  geom_point(aes(x=f_paired_mat[,1], y=f_paired_mat[,2]), size=0.5, alpha = 0.4) +
  geom_point(aes(x=f_paired_mat[1,1], y=f_paired_mat[1,2]), col = "red", pch = 4, size = 3) +
  labs(x = bquote("F"[x]^pi), y = bquote("F"[z]^pi)) +
  theme(strip.background = element_rect(fill=NA),
        text = element_text(size = 8),
        strip.text = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(fill = "transparent", color = 'black', size=0.5),
        legend.position = "Bottom",
        axis.title.x = element_text(size=8, colour='black'),
        axis.title.y = element_text(size=8, colour='black'),
        axis.text.x = element_text(size=8, colour='black'),
        axis.text.y = element_text(size=8, colour='black'),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth=0.25, colour = 'black')
  )

design <- "
  12
"
p2 + p1 + plot_layout(design = design) + plot_annotation(tag_levels = "A")

ggsave('figures/supp-Fval_rev2.pdf', width=5, height=2.6, units='in')