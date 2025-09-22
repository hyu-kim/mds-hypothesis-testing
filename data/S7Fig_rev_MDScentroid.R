## computes size and centroid distances between MDS clusters
library(ggplot2)
library('dplyr')
library(cowplot)
library(scales)


th <- theme(strip.background = element_rect(fill=NA),
            panel.background = element_rect(fill = "transparent", color = NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.background = element_blank(),
            panel.border = element_rect(fill = "transparent", color = 'black', size=0.5),
            legend.title=element_blank(), 
            legend.text=element_text(size=9),
            legend.position = "right",
            axis.title.x = element_text(size=9, colour='black'),
            axis.title.y = element_text(size=9, colour='black'),
            axis.text.x = element_text(size=9, colour='black'),
            axis.text.y = element_text(size=9, colour='black'),
            axis.line = element_blank(),
            axis.ticks = element_line(linewidth=0.25, colour = 'black'),
            plot.margin=margin(l=0.5,t=0.5,unit="cm"))


params_v <- (0:10)/10


# import and process data for sim
sim_df <- data.frame(matrix(ncol=8, nrow=0)) # dataset for sim
colnames(sim_df) <- c('N', 'hyperparameter', 'rep', 'cent_dist', 'group1_eig1', 
                      'group1_eig2', 'group2_eig1', 'group2_eig2')

for(N in c(50, 100, 200, 500)){
  y_sim <- as.matrix(read.csv(sprintf('result/ScalingStudy/sim_rev-N%g-Y.csv', N)))
  
  for(r in c(1:3)){
    for(p in params_v){
      z_embed <- read.csv(sprintf('result/ScalingStudy/sim_rev_%d/sim_rev_%d-N%g-fmds-%.2f-Z.csv',r,r,N,p))
      res1 <- cov.wt(z_embed[y_sim==0,])
      res2 <- cov.wt(z_embed[y_sim==1,])
      dist <- norm(res1$center-res2$center, type="2")
      g1e1 <- eigen(res1$cov)$values[1]
      g1e2 <- eigen(res1$cov)$values[2]
      g2e1 <- eigen(res2$cov)$values[1]
      g2e2 <- eigen(res2$cov)$values[2]
      
      sim_df <- rbind(sim_df, data.frame(N = N, hyperparameter = p, rep = r, cent_dist = dist,
                                         group1_eig1 = g1e1, group1_eig2 = g1e2,
                                         group2_eig1 = g2e1, group2_eig2 = g2e2))
    }
  }
}

sim_df_stat <- sim_df %>%
  group_by(N, hyperparameter) %>%
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
y_sim <- as.matrix(read.csv('result/Evaluation/Alga/alga_2-Y.csv'))

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
p1.1 <- ggplot(data=filter(sim_df_stat, N==50)) +
  geom_point(aes(x=hyperparameter, y=cent_dist_mean)) +
  geom_errorbar(aes(x=hyperparameter, ymin=cent_dist_mean-cent_dist_std, 
                    ymax=cent_dist_mean+cent_dist_std),
                width=0.03, color='black', size=0.25) +
  geom_line(aes(x=hyperparameter, y=cent_dist_mean), size=0.25, linetype='solid') +
  labs(x ="Hyperparameter", y = "Intergroup distance") +
  scale_x_continuous(breaks=c(0:5)/5) + 
  # scale_color_manual(values=c('black','black')) +
  th

p1.2 <- ggplot(data=filter(sim_df_stat, N==50)) +
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

p1.3 <- ggplot(data=filter(sim_df_stat, N==50)) +
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


p2.1 <- ggplot(data=filter(sim_df_stat, N==100)) +
  geom_point(aes(x=hyperparameter, y=cent_dist_mean)) +
  geom_errorbar(aes(x=hyperparameter, ymin=cent_dist_mean-cent_dist_std, 
                    ymax=cent_dist_mean+cent_dist_std),
                width=0.03, color='black', size=0.25) +
  geom_line(aes(x=hyperparameter, y=cent_dist_mean), size=0.25, linetype='solid') +
  labs(x ="Hyperparameter", y = "Intergroup distance") +
  scale_x_continuous(breaks=c(0:5)/5) + 
  # scale_color_manual(values=c('black','black')) +
  th

p2.2 <- ggplot(data=filter(sim_df_stat, N==100)) +
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

p2.3 <- ggplot(data=filter(sim_df_stat, N==100)) +
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


p3.1 <- ggplot(data=filter(sim_df_stat, N==200)) +
  geom_point(aes(x=hyperparameter, y=cent_dist_mean)) +
  geom_errorbar(aes(x=hyperparameter, ymin=cent_dist_mean-cent_dist_std, 
                    ymax=cent_dist_mean+cent_dist_std),
                width=0.03, color='black', size=0.25) +
  geom_line(aes(x=hyperparameter, y=cent_dist_mean), size=0.25, linetype='solid') +
  labs(x ="Hyperparameter", y = "Intergroup distance") +
  scale_x_continuous(breaks=c(0:5)/5) + 
  # scale_color_manual(values=c('black','black')) +
  th

p3.2 <- ggplot(data=filter(sim_df_stat, N==200)) +
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

p3.3 <- ggplot(data=filter(sim_df_stat, N==200)) +
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


p4.1 <- ggplot(data=filter(sim_df_stat, N==500)) +
  geom_point(aes(x=hyperparameter, y=cent_dist_mean)) +
  geom_errorbar(aes(x=hyperparameter, ymin=cent_dist_mean-cent_dist_std, 
                    ymax=cent_dist_mean+cent_dist_std),
                width=0.03, color='black', size=0.25) +
  geom_line(aes(x=hyperparameter, y=cent_dist_mean), size=0.25, linetype='solid') +
  labs(x ="Hyperparameter", y = "Intergroup distance") +
  scale_x_continuous(breaks=c(0:5)/5) + 
  # scale_color_manual(values=c('black','black')) +
  th

p4.2 <- ggplot(data=filter(sim_df_stat, N==500)) +
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

p4.3 <- ggplot(data=filter(sim_df_stat, N==500)) +
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


p5.1 <- ggplot(data=alg_df_stat) +
  geom_point(aes(x=hyperparameter, y=cent_dist_mean)) +
  geom_errorbar(aes(x=hyperparameter, ymin=cent_dist_mean-cent_dist_std, 
                    ymax=cent_dist_mean+cent_dist_std),
                width=0.03, color='black', size=0.25) +
  geom_line(aes(x=hyperparameter, y=cent_dist_mean), size=0.25, linetype='solid') +
  labs(x ="Hyperparameter", y = "Intergroup distance") +
  scale_x_continuous(breaks=c(0:5)/5) + 
  scale_y_continuous(labels = scientific_format(digits = 2)) +
  th

p5.2 <- ggplot(data=alg_df_stat) +
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
  scale_y_continuous(labels = scientific_format(digits = 2)) +
  th

p5.3 <- ggplot(data=alg_df_stat) +
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
  scale_y_continuous(labels = scientific_format(digits = 2)) +
  th

plot_grid(p1.1, p1.2, p1.3, p2.1, p2.2, p2.3, p3.1, p3.2, p3.3, 
          p4.1, p4.2, p4.3, p5.1, p5.2, p5.3,
          labels='AUTO', ncol=3)

ggsave('figures/Fig_S7_rev.pdf', width=6.5, height=8, units='in')



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