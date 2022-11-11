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
f_mat = matrix(0, nrow=1000, ncol=4)
for (iter in 1:1000){
  random_index <- sample(1:36, 18, replace=F)
  y <- rep(1,36)
  y[random_index] = 2  # crate random labels for ps dataset
  f_mat[iter,1] = pseudo_F(mat = distmat1, trt = y)$pseudoF
  f_mat[iter,2] = pseudo_F(mat = ordu1$vectors[,1:2], trt = y)$pseudoF
  f_mat[iter,3] = pseudo_F(mat = distmat2, trt = y)$pseudoF
  f_mat[iter,4] = pseudo_F(mat = ordu2$vectors[,1:2], trt = y)$pseudoF
}

# plot
par(mfrow = c(1,2))
plot(f_mat[,c(1,2)], main = "Site 1", 'xlab' = 'Original', 'ylab'='2D')
plot(f_mat[,c(3,4)], main = "Site 2", 'xlab' = 'Original', 'ylab'='2D')