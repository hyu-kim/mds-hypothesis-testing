library(vegan)
library(tsne)
library(uwot)

# simulated dataset
r <- 3
simrev_dist <- as.matrix(vegdist(t(readRDS(sprintf(
  'result/ScalingStudy/sim_rev_%g/sim_rev_%g-N200-data.Rds', r, r))), method="bray"))
simrev_y <- as.matrix(read.csv(sprintf('result/ScalingStudy/sim_rev-N200-Y.csv')))

alga_dist <- readRDS('Data/Alga/alga-dist.rds')


### MDS or PCoA
# simulated
z0 <- cmdscale(simrev_dist, k = 2)
write.csv(z0, sprintf('result/Evaluation/sim_rev_%g-PCoA-Z.csv', r), row.names=FALSE)



### Isomap
# simulated
for(k in c(5,10,20,50)){
  res <- isomap(dist = simrev_dist, ndim = 2, k = k)$points
  write.csv(res, sprintf('result/Evaluation/sim_rev_%g-iso-%02d-Z.csv',r, k), row.names=FALSE)
}

# algal microbiome
for(k in c(5,7,10)){
  res <- isomap(dist = data_dist, ndim = 2, k = k)$points
  write.csv(res, sprintf('Data/Alga/Isomap/alga_2-iso-%02d-Z.csv',k), row.names=FALSE)
}



### t-SNE
# simulated dataset
for(perp in c(5,10,20,50)){
  res <- tsne(X = simrev_dist, perplexity = perp)
  write.csv(res, sprintf('result/Evaluation/sim_rev_%g-tsne-%02d-Z.csv', r, perp), row.names=FALSE)
}

# algal microbiome
data_dist <- readRDS('Data/Alga/alga-dist.rds')

for(perp in c(5,7,10)){
  res <- tsne(X = dist_mat, perplexity = perp)
  write.csv(res, sprintf('Data/Alga/t-SNE/alga_2-tsne-%02d-Z.csv', perp), row.names=FALSE)
}



### UMAP
## n_neighbors: The size of local neighborhood (in terms of number of 
## neighboring sample points) used for manifold approximation. 
## Larger values result in more global views of the manifold, 
## while smaller values result in more local data being preserved. 

## y	: Optional target data for supervised dimension reduction.

# simulated dataset
for(k in c(5,10,20,50)){
  res <- umap(simrev_dist, n_neighbors = k, y = factor(simrev_y))
  write.csv(res, sprintf('result/Evaluation/sim_rev_%g-umap_s-%02d-Z.csv',r,k), row.names=FALSE)
  res <- umap(simrev_dist, n_neighbors = k)
  write.csv(res, sprintf('result/Evaluation/sim_rev_%g-umap_u-%02d-Z.csv',r,k), row.names=FALSE)
}

# algal microbiome
data_dist <- readRDS('Data/Alga/alga-dist.rds')

for(k in c(5,7,10)){
  res <- umap(data_dist, n_neighbors = k, y = factor(y))
  write.csv(res, sprintf('Data/Alga/UMAP-S/alga-umap_s-%02d-Z.csv',k), row.names=FALSE)
  res <- umap(data_dist, n_neighbors = k)
  write.csv(res, sprintf('Data/Alga/UMAP-U/alga-umap_u-%02d-Z.csv',k), row.names=FALSE)
}