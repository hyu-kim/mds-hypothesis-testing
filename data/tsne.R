library(tsne)

# simulated dataset 1
data_df <- as.matrix(read.csv(sprintf('result/Evaluation/sim_1-data.csv')))
dist_mat <- dist(data_df[,1:4])

for(perp in c(5,10,20,30)){
  res <- tsne(X = dist_mat, perplexity = perp)
  write.csv(res, sprintf('result/Evaluation/sim_1-tsne-%02d-Z.csv', perp), row.names=FALSE)
}


# algal microbiome
data_dist <- readRDS('result/Evaluation/alga_2-dist.rds')

for(perp in c(5,7,10)){
  res <- tsne(X = dist_mat, perplexity = perp)
  write.csv(res, sprintf('result/Evaluation/alga_2-iso-%02d-Z.csv', perp), row.names=FALSE)
}


## Cirrhosis and t2d
for(perp in c(5:10)){
  cirr_res$tsne[[perp]] <- tsne(X = phyl_unifrac_cirrhosis, perplexity = perp)
  t2d_res$tsne[[perp]] <- tsne(X = phyl_unifrac_t2d, perplexity = perp)
}

### F and p value
# tsne1_res <- tsne2_res <- data.frame(
#   perplexity = c(5:10),
#   pseudo_F = rep(0, 6),
#   p_val = rep(0, 6)
# )
# 
# set.seed(100)
# for(i in 1:6){
#   ttt <- get_p(mat = tsne1[[i+4]], trt = y1)
#   tsne1_res$pseudo_F[i] <- ttt$ratio
#   tsne1_res$p_val[i] <- ttt$p
#   ttt <- get_p(mat = tsne2[[i+4]], trt = y2)
#   tsne2_res$pseudo_F[i] <- ttt$ratio
#   tsne2_res$p_val[i] <- ttt$p
# }
# 
# #### Plots
# pdf("tsne.pdf", width = 9, height = 6)
# par(mfrow = c(2,3))
# for(perp in c(5:10)){
#   plot(tsne1[[perp]], col = y1, main = paste("S1, perp=", perp),
#        xlab = "config1", ylab = "config2")
#   mtext(paste("F=", round(tsne1_res$pseudo_F[perp-4], 2),
#               ", p=", tsne1_res$p_val[perp-4]), side=3)
# }
# for(perp in c(5:10)){
#   plot(tsne2[[perp]], col = y2, main = paste("S2, perp=", perp),
#        xlab = "config1", ylab = "config2")
#   mtext(paste("F=", round(tsne2_res$pseudo_F[perp-4], 2),
#               ", p=", tsne2_res$p_val[perp-4]), side=3)
# }
# dev.off()
# 
# 
# #### Shepard Plot
# ## Original distance
# pdf("tsne_shepard_orig.pdf", width = 9, height = 6)
# par(mfrow = c(2,3))
# for(perp in 5:10){
#   plot(dist1, tsne1_dist[[perp]],
#        xlab = "original distance", ylab = "configuration distance",
#        main = paste("S1, Perp=", perp, by = ""))
#   mtext(paste(
#     "Stress1 = ", 
#     round(sqrt(sum((dist1 - tsne1_dist[[perp]])^2) /
#                  sum((tsne1_dist[[perp]])^2)), 4)), side=3)
# }
# for(perp in 5:10){
#   plot(dist2, tsne2_dist[[perp]],
#        xlab = "original distance", ylab = "configuration distance",
#        main = paste("S2, Perp=", perp, by = ""))
#   mtext(paste(
#     "Stress1 = ", 
#     round(sqrt(sum((dist2 - tsne2_dist[[perp]])^2) /
#                  sum((tsne2_dist[[perp]])^2)), 4)), side=3)
# }
# dev.off()
# 
# ## Scaled distance
# pdf("tsne_shepard_scaled.pdf", width = 9, height = 6)
# par(mfrow = c(2,3))
# for(perp in 5:10){
#   plot(dist1/max(dist1), tsne1_dist[[perp]]/max(tsne1_dist[[perp]]),
#        xlab = "original distance", ylab = "configuration distance",
#        main = paste("S1,Perp=", perp, by = ""))
#   mtext(paste(
#     "Stress1 = ", 
#     round(sqrt(sum((dist1/max(dist1) - tsne1_dist[[perp]]/max(tsne1_dist[[perp]]))^2) /
#                  sum((tsne1_dist[[perp]]/max(tsne1_dist[[perp]]))^2)), 4)), side=3)
# }
# 
# par(mfrow = c(2,3))
# for(perp in 5:10){
#   plot(dist2/max(dist2), tsne2_dist[[perp]]/max(tsne2_dist[[perp]]),
#        xlab = "original distance", ylab = "configuration distance",
#        main = paste("S2,Perp=", perp, by = ""))
#   mtext(paste(
#     "Stress1 = ", 
#     round(sqrt(sum((dist2/max(dist2) - tsne2_dist[[perp]]/max(tsne2_dist[[perp]]))^2) /
#                  sum((tsne2_dist[[perp]]/max(tsne2_dist[[perp]]))^2)), 4)), side=3)
# }
# 
# dev.off()