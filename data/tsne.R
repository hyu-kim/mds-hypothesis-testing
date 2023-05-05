library(tsne)
dist1

tsne1 <- list()
for(perp in c(5:10)){
  tsne1[[perp]] <- tsne(X = dist1, perplexity = perp)
}

tsne1_dist <- list()
for(perp in c(5:10)){
  tsne1_dist[[perp]] <- dist(tsne1[[perp]])
}


tsne2 <- list()
for(perp in c(5:10)){
  tsne2[[perp]] <- tsne(X = dist2, perplexity = perp)
}

tsne2_dist <- list()
for(perp in c(5:10)){
  tsne2_dist[[perp]] <- dist(tsne2[[perp]])
}


pdf("tsne.pdf", width = 10, height = 6)
par(mfrow = c(2,3), mar = c(3,3,2,1))
for(perp in c(5:10)){
  plot(tsne1[[perp]], col = y1, main = perp)
}

par(mfrow = c(2,3), mar = c(3,3,2,1))
for(perp in c(5:10)){
  plot(tsne2[[perp]], col = y2, main = perp)
}
dev.off()


#### Shepard Plot
# Some definitions switches x and y axis
pdf("tsne_shepard.pdf", width = 10, height = 6)
## Original distance
par(mfrow = c(2,3))
for(perp in 5:10){
  plot(dist1, tsne1_dist[[perp]],
       xlab = "original distance", ylab = "configuration distance",
       main = paste("Site1,Perp=", perp, by = ""))
  mtext(paste(
    "Stress1 = ", 
    round(sqrt(sum((dist1 - tsne1_dist[[perp]])^2) /
                 sum((tsne1_dist[[perp]])^2)), 4)), side=3)
}


## Scaled distance
par(mfrow = c(2,3))
for(perp in 5:10){
  plot(dist1/max(dist1), tsne1_dist[[perp]]/max(tsne1_dist[[perp]]),
       xlab = "original distance", ylab = "configuration distance",
       main = paste("Site1,Perp=", perp, by = ""))
  mtext(paste(
    "Stress1 = ", 
    round(sqrt(sum((dist1/max(dist1) - tsne1_dist[[perp]]/max(tsne1_dist[[perp]]))^2) /
                 sum((tsne1_dist[[perp]]/max(tsne1_dist[[perp]]))^2)), 4)), side=3)
}

## Original distance
par(mfrow = c(2,3))
for(perp in 5:10){
  plot(dist2, tsne2_dist[[perp]],
       xlab = "original distance", ylab = "configuration distance",
       main = paste("Site2,Perp=", perp, by = ""))
  mtext(paste(
    "Stress1 = ", 
    round(sqrt(sum((dist2 - tsne2_dist[[perp]])^2) /
                 sum((tsne2_dist[[perp]])^2)), 4)), side=3)
}


## Scaled distance
par(mfrow = c(2,3))
for(perp in 5:10){
  plot(dist2/max(dist1), tsne2_dist[[perp]]/max(tsne2_dist[[perp]]),
       xlab = "original distance", ylab = "configuration distance",
       main = paste("Site1,Perp=", perp, by = ""))
  mtext(paste(
    "Stress1 = ", 
    round(sqrt(sum((dist2/max(dist2) - tsne2_dist[[perp]]/max(tsne2_dist[[perp]]))^2) /
                 sum((tsne2_dist[[perp]]/max(tsne2_dist[[perp]]))^2)), 4)), side=3)
}

dev.off()
## Stress-1: Scaled Stress
## For its definition, see, for example,
## Andreas Buja, Deborah F Swayne, Michael L Littman, Nathaniel Dean, Heike Hofmann & Lisha Chen (2008) 
## Data Visualization With Multidimensional Scaling, 
## Journal of Computational and Graphical Statistics, 17:2, 444-472
## (it is in our shared google drive)
