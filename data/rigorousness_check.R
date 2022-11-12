# Randomly assigns a binary label to a given phyloseq dataset theb calculates
# pseudo-F under original and MDS ordination respectively, which is repeated
# for a number of times
 
library("phyloseq")
source("permanova_with_config.R")

ps = readRDS("community_phyloseq.Rds")
ps1 = subset_samples(ps, Site=="1")
ps2 = subset_samples(ps, Site=="2")
dist1 = phyloseq::distance(ps1, method = "unifrac", weighted = T)
dist2 = phyloseq::distance(ps2, method = "unifrac", weighted = T)
distmat1 = as.matrix(dist1)
distmat2 = as.matrix(dist2)
ordu1 = ordinate(ps1, "PCoA", distance = "unifrac", weighted=TRUE)
ordu2 = ordinate(ps2, "PCoA", distance = "unifrac", weighted=TRUE)


# iterate
f_mat = matrix(0, nrow=1000, ncol=8)
random_index_set <- matrix(0, nrow=1000, ncol=18)
for (iter in 1:1000){
  random_index <- sample(1:36, 18, replace=F)
  random_index_set[iter,] <- random_index
  y <- rep(1,36)
  y[random_index] = 2  # crate random labels for ps dataset
  f_mat[iter,1] = pseudo_F(mat = distmat1, trt = y)$pseudoF
  f_mat[iter,2] = pseudo_F(mat = ordu1$vectors[,1:2], trt = y)$pseudoF
  f_mat[iter,3] = pseudo_F(mat = distmat2, trt = y)$pseudoF
  f_mat[iter,4] = pseudo_F(mat = ordu2$vectors[,1:2], trt = y)$pseudoF
}
for (c in 1:4){
  f_sorted = sort(f_mat[,c], decreasing = TRUE)
  for (iter in 1:1000){
    f_mat[iter, c+4] = which(f_mat[iter, c] == f_sorted)
  }
}
f_mat[,5:8] <- -log10(f_mat[,5:8]/1000)


# pick and save abnormal cases
case_logi <- (f_mat[,5]>2)|(f_mat[,6]>2)|(f_mat[,7]>2)|(f_mat[,8]>2)
f_mat_sub <- f_mat[case_logi,]
extremes <- matrix(0, nrow=2, ncol=2)
rownames(extremes) <- c('max', 'min')
colnames(extremes) <- c('Site 1', 'Site 2')
for (c in 1:2){
  p_ratio <- f_mat_sub[,c+4]/f_mat_sub[,c+5]
  extremes[,c] <- which(p_ratio==max(p_ratio) | p_ratio==min(p_ratio))
}
y_crit_11 <- y_crit_12 <- y_crit_21 <- y_crit_22 <- rep(1,36)
y_crit_11[random_index_set[extremes[1,1],]] <- 2
y_crit_12[random_index_set[extremes[2,1],]] <- 2
y_crit_21[random_index_set[extremes[1,2],]] <- 2
y_crit_22[random_index_set[extremes[2,2],]] <- 2


# plot F statistics
par(mfrow = c(1,2))
plot(f_mat[,c(1,2)], main = "Site 1", 'xlab' = 'Original', 'ylab'='2D')
plot(f_mat[,c(3,4)], main = "Site 2", 'xlab' = 'Original', 'ylab'='2D')

# plot p-values
par(mfrow = c(1,2))
plot(f_mat[,c(5,6)], main = "Site 1", 'xlab' = 'Original', 'ylab'='2D')
plot(f_mat[,c(7,8)], main = "Site 2", 'xlab' = 'Original', 'ylab'='2D')

# plot MDS
par(mfrow = c(1,2))
zmds2 <- plot_ordination(ps2, ordu2, axes = 1:10)$data[, 1:2]
plot(zmds1, col = y_crit_21, 
     main = paste("log_p0=", round(f_mat[extremes[1,2],7], 3),
                  "log_p_mds=", round(f_mat[extremes[1,2],8], 3)))
plot(zmds1, col = y_crit_22, 
     main = paste("log_p0=", round(f_mat[extremes[2,2],7], 3),
                  "log_p_mds=", round(f_mat[extremes[2,2],8], 3)))