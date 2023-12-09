library(uwot)
## n_neighbors: The size of local neighborhood (in terms of number of 
## neighboring sample points) used for manifold approximation. 
## Larger values result in more global views of the manifold, 
## while smaller values result in more local data being preserved. 

## y	: Optional target data for supervised dimension reduction.

n_nbr <- c(5, 10, 20, 30)
umap1 <- vector("list", length = 8)

for(i in 1:4){
  umap1[[i]] <- umap(dist1, n_neighbors = n_nbr[i], y = factor(y1))
  umap1[[i+4]] <- umap(dist1, n_neighbors = n_nbr[i])
}

umap2 <- vector("list", length = 8)
for(i in 1:4){
  umap2[[i]] <- umap(dist2, n_neighbors = n_nbr[i], y = factor(y2))
  umap2[[i+4]] <- umap(dist2, n_neighbors = n_nbr[i])
}

### F and p value
umap1_res <- umap2_res <- data.frame(
  Supervised = c(rep(1, 4), rep(0, 4)),
  n_neighbors = c(n_nbr, n_nbr),
  pseudo_F = rep(0, 8),
  p_val = rep(0, 8)
)

set.seed(100)
for(i in 1:8){
  ttt <- get_p(mat = umap1[[i]], trt = y1)
  umap1_res$pseudo_F[i] <- ttt$ratio
  umap1_res$p_val[i] <- ttt$p
  ttt <- get_p(mat = umap2[[i]], trt = y2)
  umap2_res$pseudo_F[i] <- ttt$ratio
  umap2_res$p_val[i] <- ttt$p
}

### distance matrix
umap1_dist <- umap2_dist <- vector("list", length = 8)
for(i in c(1:8)){
  umap1_dist[[i]] <- dist(umap1[[i]])
  umap2_dist[[i]] <- dist(umap2[[i]])
}

############# Plots
pdf("umap.pdf", width = 12, height = 6)
par(mfrow = c(2,4))
for(i in 1:4){
  plot(umap1[[i]], col = y1, main = paste("S1, n=", n_nbr[i], ", Sup"),
       xlab = "config1", ylab = "config2")
  mtext(paste("F=", round(umap1_res$pseudo_F[i], 2),
              ", p=", umap1_res$p_val[i]), side=3)
}
for(i in 1:4){
  plot(umap1[[i+4]], col = y1, main = paste("S1, n=", n_nbr[i], ", UnSup"),
       xlab = "config1", ylab = "config2")
  mtext(paste("F=", round(umap1_res$pseudo_F[i+4], 2),
              ", p=", umap1_res$p_val[i+4]), side=3)
}
for(i in 1:4){
  plot(umap2[[i]], col = y2, main = paste("S2, n=", n_nbr[i], ", Sup"),
       xlab = "config1", ylab = "config2")
  mtext(paste("F=", round(umap2_res$pseudo_F[i], 2),
              ", p=", umap2_res$p_val[i]), side=3)
}
for(i in 1:4){
  plot(umap2[[i+4]], col = y2, main = paste("S2, n=", n_nbr[i], ", UnSup"),
       xlab = "config1", ylab = "config2")
  mtext(paste("F=", round(umap2_res$pseudo_F[i+4], 2),
              ", p=", umap2_res$p_val[i+4]), side=3)
}
dev.off()


#### Shepard Plot
# original
pdf("umap_shepard_orig.pdf", width = 12, height = 6)
par(mfrow = c(2,4))
for(i in 1:4){
  plot(dist1, umap1_dist[[i]], 
       xlab = "original distance", ylab = "configuration distance",
       main = paste("S1, n=", n_nbr[i], ", Sup"))
  mtext(paste(
    "Stress1 = ", 
    round(sqrt(sum((dist1 - umap1_dist[[i]])^2) /
                 sum((umap1_dist[[i]])^2)), 4)), side=3)
}
for(i in 1:4){
  plot(dist1, umap1_dist[[i+4]], 
       xlab = "original distance", ylab = "configuration distance",
       main = paste("S1, n=", n_nbr[i], ", UnSup"))
  mtext(paste(
    "Stress1 = ", 
    round(sqrt(sum((dist1 - umap1_dist[[i+4]])^2) /
                 sum((umap1_dist[[i+4]])^2)), 4)), side=3)
}
for(i in 1:4){
  plot(dist2, umap2_dist[[i]], 
       xlab = "original distance", ylab = "configuration distance",
       main = paste("S2, n=", n_nbr[i], ", Sup"))
  mtext(paste(
    "Stress1 = ", 
    round(sqrt(sum((dist2 - umap2_dist[[i]])^2) /
                 sum((umap2_dist[[i]])^2)), 4)), side=3)
}
for(i in 1:4){
  plot(dist2, umap2_dist[[i+4]], 
       xlab = "original distance", ylab = "configuration distance",
       main = paste("S2, n=", n_nbr[i], ", UnSup"))
  mtext(paste(
    "Stress1 = ", 
    round(sqrt(sum((dist2 - umap2_dist[[i+4]])^2) /
                 sum((umap2_dist[[i+4]])^2)), 4)), side=3)
}
dev.off()

## scaled
pdf("umap_shepard_scaled.pdf", width = 12, height = 6)
par(mfrow = c(2,4))
for(i in 1:4){
  plot(dist1/max(dist1), umap1_dist[[i]]/max(umap1_dist[[i]]), 
       xlab = "original distance", ylab = "configuration distance",
       main = paste("S1, n=", n_nbr[i], ", Sup"))
  mtext(paste(
    "Stress1 = ", 
    round(sqrt(sum((dist1/max(dist1) - umap1_dist[[i]]/max(umap1_dist[[i]]))^2) /
                 sum((umap1_dist[[i]]/max(umap1_dist[[i]]))^2)), 4)), side=3)
}
for(i in 1:4){
  plot(dist1/max(dist1), umap1_dist[[i+4]]/max(umap1_dist[[i+4]]), 
       xlab = "original distance", ylab = "configuration distance",
       main = paste("S1, n=", n_nbr[i], ", UnSup"))
  mtext(paste(
    "Stress1 = ", 
    round(sqrt(sum((dist1/max(dist1) - umap1_dist[[i+4]]/max(umap1_dist[[i+4]]))^2) /
                 sum((umap1_dist[[i+4]]/max(umap1_dist[[i+4]]))^2)), 4)), side=3)
}
for(i in 1:4){
  plot(dist2/max(dist2), umap2_dist[[i]]/max(umap2_dist[[i]]), 
       xlab = "original distance", ylab = "configuration distance",
       main = paste("S2, n=", n_nbr[i], ", Sup"))
  mtext(paste(
    "Stress1 = ", 
    round(sqrt(sum((dist2/max(dist2) - umap2_dist[[i]]/max(umap2_dist[[i]]))^2) /
                 sum((umap2_dist[[i]]/max(umap2_dist[[i]]))^2)), 4)), side=3)
}
for(i in 1:4){
  plot(dist2/max(dist2), umap2_dist[[i+4]]/max(umap2_dist[[i+4]]), 
       xlab = "original distance", ylab = "configuration distance",
       main = paste("S2, n=", n_nbr[i], ", UnSup"))
  mtext(paste(
    "Stress1 = ", 
    round(sqrt(sum((dist2/max(dist2) - umap2_dist[[i+4]]/max(umap2_dist[[i+4]]))^2) /
                 sum((umap2_dist[[i+4]]/max(umap2_dist[[i+4]]))^2)), 4)), side=3)
}
dev.off()
