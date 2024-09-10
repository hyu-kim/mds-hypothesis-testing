## Evaluation
### RF-PHATE paper
library("FastKNN")
knn_classifier <- function(dm, y, k){
  n <- nrow(dm)
  predictions <- rep(0, n)
  for (i in 1:n) {
    indices <- k.nearest.neighbors(i, dm, k = k)
    vec_mas_cer <- table(y[indices])
    predictions[i] <- sample(names(which(vec_mas_cer == max(vec_mas_cer))), 
                             size = 1)
  }
  tab <- table(predictions, y)
  acc <- (tab[1,1] + tab[2,2])/n
  return(acc)
}
knn_classifier(distmat1, y1, 5)

knn_regression <- function(dm, emb, k){
  #emb: n * 2 matrix
  n <- nrow(dm)
  predictions <- matrix(0, nrow = n, ncol = ncol(emb))
  for (i in 1:n) {
    indices <- k.nearest.neighbors(i, dm, k = k)
    km <- apply(emb[indices,], 2, mean)
    predictions[i, ] <- km
  }
  mse <- mean(apply(emb - predictions, 1, function(x){sum(x^2)}))
  return(mse = mse)
}
knn_regression(distmat1, res1$lambda0.5$z, 5)


rp_eval <- function(nperm, dm, y, emb, k){
  n <- length(y)
  rho1 <- rep(knn_classifier(dm, y, k), nperm) #ACC: higher better
  rho2 <- rep(knn_regression(dm, emb, k), nperm) #MSE: lower better
  for(i in 1:nperm){
    set.seed(i)
    ind_perm <- sample(1:n, size = n, replace = F)
    y_perm <- y[ind_perm]
    emb_perm <- emb[ind_perm,]
    rho1[i] <- rho1[i] - knn_classifier(dm, y_perm, k)
    rho2[i] <- knn_regression(dm, emb_perm, k) - rho2[i]
  }
  return(cor(rho1, rho2))
}



rp_eval_n_sim <- rep(0, 100)
set.seed(100)
for(i in 1:100){
  tt <- rp_eval(200, sim_data$distmat, sim_data$data$Y, sim_res$proposed$lambda0.5$z, 
                5)
  rp_eval_n_sim[i] <- tt$diff_corr
}
rp_eval_n_sim <- max(rp_eval_n_sim) - rp_eval_n_sim
summary(rp_eval_n_sim)


### Isomap
library(vegan)
isomap1 <- list(
  isomap(dist = dist1, ndim = 2, k = 5),
  isomap(dist = dist1, ndim = 2, k = 7),
  isomap(dist = dist1, ndim = 2, k = 10))
plot(isomap1[[3]]$points)

isomap2 <- list(
  isomap(dist = dist2, ndim = 2, k = 5),
  isomap(dist = dist2, ndim = 2, k = 7),
  isomap(dist = dist2, ndim = 2, k = 10))
plot(isomap2[[3]]$points)

sim_res$isomap <- list(
  isomap(dist = sim_data$dist, ndim = 2, k = 5),
  isomap(dist = sim_data$dist, ndim = 2, k = 10),
  isomap(dist = sim_data$dist, ndim = 2, k = 20),
  isomap(dist = sim_data$dist, ndim = 2, k = 30)
  )

## tsne
library(tsne)
sim_res$tsne <- list(
  tsne(X = sim_data$dist, perplexity = 5),
  tsne(X = sim_data$dist, perplexity = 10),
  tsne(X = sim_data$dist, perplexity = 20),
  tsne(X = sim_data$dist, perplexity = 30)
) 


### PHATE : some python module error
# library(phateR)
# phate(data = distmat1,
#       ndim = 2,
#       knn.dist.method = "precomputed")
# data(tree.data.small)
# phate(tree.data.small$data)

### Eval
eval_res <- list()
# eval_res$mds <- rep(0, 100)
# eval_res$fmds3 <- rep(0, 100)
# eval_res$fmds5 <- rep(0, 100)
# eval_res$fmds7 <- rep(0, 100)
# 
# eval_res$smds <- rep(0, 100)
# # UMAP
# eval_res$umap_uns5 <- rep(0, 100)
# eval_res$umap_uns10 <- rep(0, 100)
# eval_res$umap_uns20 <- rep(0, 100)
# eval_res$umap_uns30 <- rep(0, 100)
# eval_res$umap_sup5 <- rep(0, 100)
# eval_res$umap_sup10 <- rep(0, 100)
# eval_res$umap_sup20 <- rep(0, 100)
# eval_res$umap_sup30 <- rep(0, 100)
# # tsne
# eval_res$tsne5 <- rep(0, 100)
# eval_res$tsne10 <- rep(0, 100)
# eval_res$tsne20 <- rep(0, 100)
# eval_res$tsne30 <- rep(0, 100)
# # Isomap
# eval_res$isomap5 <- rep(0, 100)
# eval_res$isomap10 <- rep(0, 100)
# eval_res$isomap20 <- rep(0, 100)
# eval_res$isomap30 <- rep(0, 100)
# 
# sim_res$proposed$lambda0.7$z
# sim_res2$proposed$lambda0.7$z
# 
# for(i in 1:100){
#   # mds
#   set.seed(i)
#   eval_res$mds[i] <- rp_eval(100, sim_data$distmat, sim_data$data$Y, 
#                 sim_res$mds, 10)$diff_corr
#   # fmds
#   set.seed(i)
#   eval_res$fmds3[i] <- rp_eval(100, sim_data$distmat, sim_data$data$Y, 
#                                sim_res$proposed$lambda0.3$z, 10)$diff_corr
#   set.seed(i)
#   eval_res$fmds5[i] <- rp_eval(100, sim_data$distmat, sim_data$data$Y, 
#                                sim_res$proposed$lambda0.5$z, 10)$diff_corr
#   set.seed(i)
#   eval_res$fmds7[i] <- rp_eval(100, sim_data$distmat, sim_data$data$Y, 
#                                sim_res$proposed$lambda0.7$z, 10)$diff_corr
#   # smds
#   
#   
#   
#   # UMAP Unsupervised
#   set.seed(i)
#   eval_res$umap_uns5[i] <- rp_eval(100, sim_data$distmat, sim_data$data$Y, 
#                                    sim_res$umap[[1]], 10)$diff_corr
#   set.seed(i)
#   eval_res$umap_uns10[i] <- rp_eval(100, sim_data$distmat, sim_data$data$Y, 
#                                    sim_res$umap[[2]], 10)$diff_corr
#   set.seed(i)
#   eval_res$umap_uns20[i] <- rp_eval(100, sim_data$distmat, sim_data$data$Y, 
#                                    sim_res$umap[[3]], 10)$diff_corr
#   set.seed(i)
#   eval_res$umap_uns30[i] <- rp_eval(100, sim_data$distmat, sim_data$data$Y, 
#                                    sim_res$umap[[4]], 10)$diff_corr
#   # UMAP Supervised
#   set.seed(i)
#   eval_res$umap_sup5[i] <- rp_eval(100, sim_data$distmat, sim_data$data$Y, 
#                                    sim_res$umap[[5]], 10)$diff_corr
#   set.seed(i)
#   eval_res$umap_sup10[i] <- rp_eval(100, sim_data$distmat, sim_data$data$Y, 
#                                     sim_res$umap[[6]], 10)$diff_corr
#   set.seed(i)
#   eval_res$umap_sup20[i] <- rp_eval(100, sim_data$distmat, sim_data$data$Y, 
#                                     sim_res$umap[[7]], 10)$diff_corr
#   set.seed(i)
#   eval_res$umap_sup30[i] <- rp_eval(100, sim_data$distmat, sim_data$data$Y, 
#                                     sim_res$umap[[8]], 10)$diff_corr
#   # tsne
#   set.seed(i)
#   eval_res$tsne5[i] <- rp_eval(100, sim_data$distmat, sim_data$data$Y, 
#                                sim_res$tsne[[1]], 10)$diff_corr
#   set.seed(i)
#   eval_res$tsne10[i] <- rp_eval(100, sim_data$distmat, sim_data$data$Y, 
#                                sim_res$tsne[[2]], 10)$diff_corr
#   set.seed(i)
#   eval_res$tsne20[i] <- rp_eval(100, sim_data$distmat, sim_data$data$Y, 
#                                sim_res$tsne[[3]], 10)$diff_corr
#   set.seed(i)
#   eval_res$tsne30[i] <- rp_eval(100, sim_data$distmat, sim_data$data$Y, 
#                                sim_res$tsne[[4]], 10)$diff_corr
#   # Isomap
#   set.seed(i)
#   eval_res$isomap5[i] <- rp_eval(100, sim_data$distmat, sim_data$data$Y, 
#                                  sim_res$isomap[[1]]$points, 10)$diff_corr
#   set.seed(i)
#   eval_res$isomap10[i] <- rp_eval(100, sim_data$distmat, sim_data$data$Y, 
#                                   sim_res$isomap[[2]]$points, 10)$diff_corr
#   set.seed(i)
#   eval_res$isomap20[i] <- rp_eval(100, sim_data$distmat, sim_data$data$Y, 
#                                   sim_res$isomap[[3]]$points, 10)$diff_corr
#   set.seed(i)
#   eval_res$isomap30[i] <- rp_eval(100, sim_data$distmat, sim_data$data$Y, 
#                                   sim_res$isomap[[4]]$points, 10)$diff_corr
# }



lapply(eval_res, summary)

lapply(eval_res, mean)

for(i in 1:100){
  # fmds
  set.seed(i)
  eval_res$fmds3[i] <- rp_eval(100, sim_data$distmat, sim_data$data$Y, 
                               sim_res2$proposed$lambda0.3$z, 10)$diff_corr
  set.seed(i)
  eval_res$fmds5[i] <- rp_eval(100, sim_data$distmat, sim_data$data$Y, 
                               sim_res2$proposed$lambda0.5$z, 10)$diff_corr
  set.seed(i)
  eval_res$fmds7[i] <- rp_eval(100, sim_data$distmat, sim_data$data$Y, 
                               sim_res$proposed$lambda0.7$z, 10)$diff_corr
}
# 
# eval_res$fmds3 <- max(eval_res$fmds3) - eval_res$fmds3
# summary(rp_eval_n)

######### eval to site1
eval_res$site1 <- lapply(
  site1_embed[1:21],
  function(x){
    rp_eval(100, distmat1, y1, x, 10)
  }
)
## without nn

eval_res$site2 <- lapply(
  site2_embed[1:21],
  function(x){
    rp_eval(100, distmat2, y2, x, 10)
  }
)

eval_res$sim <- lapply(
  sim_embed,
  function(x){
    rp_eval(100, sim_data$distmat, sim_data$data$Y, x, 10)
  }
)


cbind(eval_res$site1, eval_res$site2)
t(data.frame(eval_res$sim))
