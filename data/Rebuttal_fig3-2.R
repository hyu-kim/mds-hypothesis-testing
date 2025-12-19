## computes F-MDS w/o convergence limit
library(vegan)
library(parallel)
source("mm.R")

run_mm_cmds <- function(params) {
  require(vegan)
  require("mm.R")
  require("permanova_with_config.R")
  r <- params$replicate
  N <- params$size
  l <- params$lambda
  
  data_mat <- readRDS('result/ScalingStudy/sim_rev_1/sim_rev_1-N100-data.Rds')
  dist_mat <- as.matrix(vegdist(t(data_mat), method="bray"))
  y_df <- as.matrix(read.csv(sprintf('result/ScalingStudy/sim_rev-N100-Y.csv')))
  
  print(sprintf("Data Read. N: %d, Replicate: %d, For lambda: %g", N, r, l))
  
  z0 <- cmdscale(dist_mat, k = 2)
  res <- mm_cmds(nit = 50, lambda = l, z0 = z0, D = dist_mat, y = y_df,
                 threshold = 0, dataset = 'sim_rev_1-N100', 
                 folder = "result/Revision2_Dec2025")
  
  # write.csv(res, sprintf('Simulated/F-MDS/Results/sim_rev_%d-N%d-fmds-%.2f-Z', r, N, l), row.names=FALSE)
  
  print(sprintf("MM Done. N: %d, Replicate: %d, For lambda: %g", N, r, l))
}

# dataframe of combinations
replicate_list <- c(1)
size_list <- c(100)
lambda_list <- c(0.2, 0.5, 0.8)

# loop and run parellel
param_combinations <- expand.grid(replicate_list, size_list, lambda_list)
colnames(param_combinations) <- c("replicate", "size", "lambda")

num_cores <- detectCores()
cl <- makeCluster(num_cores - 1) # Leave one core free for system tasks

clusterExport(cl, ls(envir = .GlobalEnv), envir = environment())

results <- parLapply(cl, split(param_combinations, seq_len(nrow(param_combinations))), run_mm_cmds)

stopCluster(cl)


## Plot results
log_all_df <- data.frame(matrix(nrow=0, ncol=6))
colnames(log_all_df) <- c('dataset', 'size', 'lambda', 'epoch', 'p_z', 'p_0')

for(r in c(1)){
  for(N in c(100)){
    for(l in c(0.2, 0.5, 0.8)){
      log_df <- read.csv(sprintf('result/Revision2_Dec2025/sim_rev_%g-N%g-fmds-%.2f-log.csv', r, N, l))
      log_all_df <- rbind(log_all_df, 
                          data.frame(dataset = sprintf('sim_rev_%g', r), 
                                     size = N, lambda = l, 
                                     epoch = log_df$epoch, 
                                     p_z = log_df$p_z, p_0 = log_df$p_0
                          ))
    }
  }
}

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

