## computes size and centroid distances between MDS clusters
library(ggplot2)
library('dplyr')
library(cowplot)


th <- theme(strip.background = element_rect(fill=NA),
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
            axis.ticks = element_line(linewidth=0.25, colour = 'black'))


params_v <- (0:10)/10


# import and process data for sim
sim_df <- data.frame(matrix(ncol=7, nrow=0)) # dataset for sim
colnames(sim_df) <- c('hyperparameter', 'rep', 'cent_dist', 'group1_eig1', 
                      'group1_eig2', 'group2_eig1', 'group2_eig2')
sim_data <- as.matrix(read.csv('result/Evaluation/sim_1-data.csv'))

for(r in c(1:3)){
  y_sim <- as.matrix(read.csv(sprintf('result/HyperparameterStudy/sim_%d-Y.csv', r)))
  for(p in params_v){
    z_embed <- read.csv(sprintf('result/HyperparameterStudy/sim_%d/sim_%d-fmds-%.2f-Z.csv',r,r,p))
    res1 <- cov.wt(z_embed[y_sim==1,])
    res2 <- cov.wt(z_embed[y_sim==2,])
    dist <- norm(res1$center-res2$center, type="2")
    g1e1 <- eigen(res1$cov)$values[1]
    g1e2 <- eigen(res1$cov)$values[2]
    g2e1 <- eigen(res2$cov)$values[1]
    g2e2 <- eigen(res2$cov)$values[2]
    
    sim_df <- rbind(sim_df, data.frame(hyperparameter = p, rep = r, cent_dist = dist,
                                       group1_eig1 = g1e1, group1_eig2 = g1e2,
                                       group2_eig1 = g2e1, group2_eig2 = g2e2))
  }
}

sim_df_stat <- sim_df %>%
  group_by(hyperparameter) %>%
  summarize(cent_dist_mean = mean(cent_dist), cent_dist_std = sd(cent_dist), 
            g1e1_av = mean(group1_eig1), g1e1_sd = sd(group1_eig1),
            g1e2_av = mean(group1_eig2), g1e2_sd = sd(group1_eig2),
            g2e1_av = mean(group2_eig1), g2e1_sd = sd(group2_eig1),
            g2e2_av = mean(group2_eig2), g2e2_sd = sd(group2_eig2),
            n = n())


# import and process data for alg
alg_df <- data.frame(matrix(ncol=6, nrow=0)) # dataset for alg
colnames(alg_df) <- c('hyperparameter', 'cent_dist', 'group1_eig1', 
                      'group1_eig2', 'group2_eig1', 'group2_eig2')
y_sim <- as.matrix(read.csv('result/Evaluation/alga_2-Y.csv'))

for(p in params_v){
  z_embed <- read.csv(sprintf('result/HyperparameterStudy/alga_2/alga_2-fmds-%.2f-Z.csv',p))
  res1 <- cov.wt(z_embed[y_sim==1,])
  res2 <- cov.wt(z_embed[y_sim==2,])
  dist <- norm(res1$center-res2$center, type="2")
  g1e1 <- eigen(res1$cov)$values[1]
  g1e2 <- eigen(res1$cov)$values[2]
  g2e1 <- eigen(res2$cov)$values[1]
  g2e2 <- eigen(res2$cov)$values[2]
  
  alg_df <- rbind(alg_df, data.frame(hyperparameter = p, cent_dist = dist,
                                     group1_eig1 = g1e1, group1_eig2 = g1e2,
                                     group2_eig1 = g2e1, group2_eig2 = g2e2))
}

alg_df_stat <- alg_df %>%
  group_by(hyperparameter) %>%
  summarize(cent_dist_mean = mean(cent_dist), cent_dist_std = sd(cent_dist), 
            g1e1_av = mean(group1_eig1), g1e1_sd = sd(group1_eig1),
            g1e2_av = mean(group1_eig2), g1e2_sd = sd(group1_eig2),
            g2e1_av = mean(group2_eig1), g2e1_sd = sd(group2_eig1),
            g2e2_av = mean(group2_eig2), g2e2_sd = sd(group2_eig2),
            n = n())


## Plot Figure
p1 <- ggplot(data=sim_df_stat) +
  geom_point(aes(x=hyperparameter, y=cent_dist_mean)) +
  geom_errorbar(aes(x=hyperparameter, ymin=cent_dist_mean-cent_dist_std, 
                    ymax=cent_dist_mean+cent_dist_std),
                width=0.03, color='black', size=0.25) +
  geom_line(aes(x=hyperparameter, y=cent_dist_mean), size=0.25, linetype='solid') +
  labs(x ="Hyperparameter", y = "Group centroid distance") +
  scale_x_continuous(breaks=c(0:5)/5) + 
  # scale_color_manual(values=c('black','black')) +
  th


p2 <- ggplot(data=sim_df_stat) +
  geom_point(aes(x=hyperparameter, y=g1e1_av), color='blue') +
  geom_errorbar(aes(x=hyperparameter, ymin=g1e1_av-g1e1_sd, ymax=g1e1_av+g1e1_sd),
                width=0.03, color='blue', size=0.25) +
  geom_line(aes(x=hyperparameter, y=g1e1_av), size=0.25, linetype='solid', color='blue') +
  
  geom_point(aes(x=hyperparameter, y=g2e1_av), color='red') +
  geom_errorbar(aes(x=hyperparameter, ymin=g2e1_av-g2e1_sd, ymax=g2e1_av+g2e1_sd),
                width=0.03, color='red', size=0.25) +
  geom_line(aes(x=hyperparameter, y=g2e1_av), size=0.25, linetype='solid', color='red') +
  labs(x ="Hyperparameter", y = "Long axis variance") +
  scale_x_continuous(breaks=c(0:5)/5) + 
  th


p3 <- ggplot(data=sim_df_stat) +
  geom_point(aes(x=hyperparameter, y=g1e2_av), color='blue') +
  geom_errorbar(aes(x=hyperparameter, ymin=g1e2_av-g1e2_sd, ymax=g1e2_av+g1e2_sd),
                width=0.03, color='blue', size=0.25) +
  geom_line(aes(x=hyperparameter, y=g1e2_av), size=0.25, linetype='solid', color='blue') +
  
  geom_point(aes(x=hyperparameter, y=g2e2_av), color='red') +
  geom_errorbar(aes(x=hyperparameter, ymin=g2e2_av-g2e2_sd, ymax=g2e2_av+g2e2_sd),
                width=0.03, color='red', size=0.25) +
  geom_line(aes(x=hyperparameter, y=g2e2_av), size=0.25, linetype='solid', color='red') +
  labs(x ="Hyperparameter", y = "Short axis variance") +
  scale_x_continuous(breaks=c(0:5)/5) + 
  th


p4 <- ggplot(data=alg_df_stat) +
  geom_point(aes(x=hyperparameter, y=cent_dist_mean)) +
  geom_errorbar(aes(x=hyperparameter, ymin=cent_dist_mean-cent_dist_std, 
                    ymax=cent_dist_mean+cent_dist_std),
                width=0.03, color='black', size=0.25) +
  geom_line(aes(x=hyperparameter, y=cent_dist_mean), size=0.25, linetype='solid') +
  labs(x ="Hyperparameter", y = "Group centroid distance") +
  scale_x_continuous(breaks=c(0:5)/5) + 
  th


p5 <- ggplot(data=alg_df_stat) +
  geom_point(aes(x=hyperparameter, y=g1e1_av), color='blue') +
  geom_errorbar(aes(x=hyperparameter, ymin=g1e1_av-g1e1_sd, ymax=g1e1_av+g1e1_sd),
                width=0.03, color='blue', size=0.25) +
  geom_line(aes(x=hyperparameter, y=g1e1_av), size=0.25, linetype='solid', color='blue') +
  
  geom_point(aes(x=hyperparameter, y=g2e1_av), color='red') +
  geom_errorbar(aes(x=hyperparameter, ymin=g2e1_av-g2e1_sd, ymax=g2e1_av+g2e1_sd),
                width=0.03, color='red', size=0.25) +
  geom_line(aes(x=hyperparameter, y=g2e1_av), size=0.25, linetype='solid', color='red') +
  labs(x ="Hyperparameter", y = "Long axis variance") +
  scale_x_continuous(breaks=c(0:5)/5) + 
  th


p6 <- ggplot(data=alg_df_stat) +
  geom_point(aes(x=hyperparameter, y=g1e2_av), color='blue') +
  geom_errorbar(aes(x=hyperparameter, ymin=g1e2_av-g1e2_sd, ymax=g1e2_av+g1e2_sd),
                width=0.03, color='blue', size=0.25) +
  geom_line(aes(x=hyperparameter, y=g1e2_av), size=0.25, linetype='solid', color='blue') +
  
  geom_point(aes(x=hyperparameter, y=g2e2_av), color='red') +
  geom_errorbar(aes(x=hyperparameter, ymin=g2e2_av-g2e2_sd, ymax=g2e2_av+g2e2_sd),
                width=0.03, color='red', size=0.25) +
  geom_line(aes(x=hyperparameter, y=g2e2_av), size=0.25, linetype='solid', color='red') +
  labs(x ="Hyperparameter", y = "Short axis variance") +
  scale_x_continuous(breaks=c(0:5)/5) + 
  th

plot_grid(p1, p2, p3, p4, p5, p6, labels='AUTO', ncol=3)

ggsave('figures/Fig_S7.pdf', width=6.5, height=4, units='in')



## Numbers to insert in main text
# increase rate of sim data for lambda<0.4 and lambda>=0.4
sim_diff_early <- sim_diff_late <- rep(0,3)
for(r in c(1:3)){
  sim_diff_early[r] <- 
    sim_df$cent_dist[sim_df$rep==r & sim_df$hyperparameter==0.3] - 
    sim_df$cent_dist[sim_df$rep==r & sim_df$hyperparameter==0]
  sim_diff_late[r] <- 
    sim_df$cent_dist[sim_df$rep==r & sim_df$hyperparameter==1] - 
    sim_df$cent_dist[sim_df$rep==r & sim_df$hyperparameter==0.4]
}
sim_diff_early <- sim_diff_early/0.3
sim_diff_late <- sim_diff_late/0.6
print(c(mean(sim_diff_early), sd(sim_diff_early), mean(sim_diff_late), sd(sim_diff_late)))

# spearman correlation of panel D
cor(sim_df_stat$hyperparameter, sim_df_stat$cent_dist_mean, method="spearman")
cor.test(sim_df_stat$hyperparameter, sim_df_stat$cent_dist_mean, method="spearman")
cor(alg_df_stat$hyperparameter, alg_df_stat$cent_dist_mean, method="spearman")
cor.test(alg_df_stat$hyperparameter, alg_df_stat$cent_dist_mean, method="spearman")

print(c(
  (alg_df_stat$g1e1_av[1] - alg_df_stat$g1e1_av[11])/alg_df_stat$g1e1_av[1],
  (alg_df_stat$g2e1_av[1] - alg_df_stat$g2e1_av[11])/alg_df_stat$g2e1_av[1],
  (alg_df_stat$g1e2_av[1] - alg_df_stat$g1e2_av[11])/alg_df_stat$g1e2_av[1],
  (alg_df_stat$g2e2_av[1] - alg_df_stat$g2e2_av[11])/alg_df_stat$g2e2_av[1]))