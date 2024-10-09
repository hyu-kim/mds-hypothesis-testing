library(ggplot2)

## functions
# Distance matrix
get_dist_mat <- function(z){
  z_dist <- dist(z)
  return(z_dist)
}

# Stress-1
get_stress <- function(x_dist, z_dist){
  res <- sqrt(sum((x_dist - z_dist)^2) / sum(z_dist^2))
  return(res)
}

# Pearson correlation
get_pearson_corr <- function(x_dist, z_dist){
  res <- cov(x_dist, z_dist) / (sd(x_dist) * sd(z_dist))
  return(res)
}


## import and process
v_lambda = (0:20)/20
df_eval <- data.frame(matrix(ncol=4, nrow=0))
colnames(df_eval) <- c('method', 'lambda', 'stress', 'pearson_corr')

for(i in 1:length(v_lambda)){
  lambda <- v_lambda[i]
  z_lambda_df <- read.csv(paste('result/simulated/iter_', 
                                ifelse(lambda>0.15, '50', '200'), '/sim-', 
                                sprintf('%.2f',lambda), '-Z.csv', sep=''))
  y_df <- read.csv(paste('result/simulated/iter_', 
                         ifelse(lambda>0.15, '50', '200'), '/sim-Y.csv', sep=''))
  x_df <- read.csv(paste('result/simulated/iter_', 
                         ifelse(lambda>0.15, '50', '200'), 
                         '/sim-data.csv', sep=''))
  x_dist <- get_dist_mat(x_df)
  z_dist <- get_dist_mat(z_lambda_df)
  stress_lambda <- get_stress(x_dist, z_dist)
  corr_lambda <- get_pearson_corr(x_dist, z_dist)
  df_eval[nrow(df_eval)+1,] <- list('FMDS', lambda, stress_lambda, corr_lambda)
}


## Panel A
ggplot() +
  geom_line(data=df_eval, aes(x=lambda, y=stress)) +
  theme(strip.background = element_rect(fill=NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(fill = "transparent", color = 'black', size=0.5),
        legend.position = "None",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=8, colour='black'),
        axis.text.y = element_text(size=8, colour='black'),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth=0.5)
  )
ggsave('result/fig2_stress_draft.pdf', width=2, height=2, units='in')

## Panel B
ggplot() +
  geom_line(data=df_eval, aes(x=lambda, y=pearson_corr))