# Randomly assigns a binary label to a given phyloseq dataset, 
# Then calculates pseudo-F under original and MDS ordination respectively, 
# which is iterated for a number of times.
 
library("phyloseq")
source("permanova_with_config.R")


# function
get_random_f_stats <- function(distmat, ordu, n_iter = 1000){
  # initialize
  random_f = matrix(0, nrow=n_iter, ncol=4)
  colnames(random_f) <- c('F_original', 'F_2D', 'p_original', 'p_2D')
  random_labels <- matrix(0, nrow=n_iter, ncol=36)
  # randomly assign labels by iteration
  print('processing iterations..')
  for (iter in 1:n_iter){
    random_index <- sample(1:36, 18, replace=F)
    y <- rep(1,36)
    y[random_index] = 2  # crate random labels for ps dataset
    random_labels[iter,] <- y
    random_f[iter,1] = pseudo_F(d = distmat, trt = y)$pseudoF
    random_f[iter,2] = pseudo_F(mat = ordu$vectors[,1:2], trt = y)$pseudoF
  }
  # compute and fill in p values
  for (c in 1:2){
    f_sorted = sort(random_f[,c], decreasing = TRUE)
    for (iter in 1:n_iter){
      random_f[iter, c+2] = which(random_f[iter, c] == f_sorted)
    }
  }
  random_f[,3:4] <- -log10(random_f[,3:4]/n_iter)
  
  return(list(f=random_f, labels=random_labels))
}


ps = readRDS("community_phyloseq.Rds")
ps1 = subset_samples(ps, Site=="1")
ps2 = subset_samples(ps, Site=="2")
dist1 = phyloseq::distance(ps1, method = "unifrac", weighted = T)
dist2 = phyloseq::distance(ps2, method = "unifrac", weighted = T)
distmat1 = as.matrix(dist1)
distmat2 = as.matrix(dist2)
ordu1 = ordinate(ps1, "PCoA", distance = "unifrac", weighted=TRUE)
ordu2 = ordinate(ps2, "PCoA", distance = "unifrac", weighted=TRUE)

f_stats1 <- get_random_f_stats(distmat=distmat1, ordu=ordu1, n_iter=2500)
f_stats2 <- get_random_f_stats(distmat=distmat2, ordu=ordu2, n_iter=2500)
random_f1 <- f_stats1$f
random_f2 <- f_stats2$f
labels1 <- f_stats1$labels
labels2 <- f_stats2$labels

# plot F statistics
par(mfrow = c(1,2))
plot(random_f1[,1:2], main = "Site 1", 'xlab' = 'pseudo F, original', 'ylab'='pseudo F, 2D')
plot(random_f2[,1:2], main = "Site 2", 'xlab' = 'pseudo F, original', 'ylab'='pseudo F, 2D')

# plot p-values
par(mfrow = c(1,2))
plot(random_f1[,3:4], main = "Site 1", 'xlab' = '-log P, original', 'ylab'='-log P, 2D')
plot(random_f2[,3:4], main = "Site 2", 'xlab' = '-log P, original', 'ylab'='-log P, 2D')