source('fig_util.R')
# import all data by Soobin
load('/Users/hk/Downloads/all_results.RData')

## SIMULATION
# export data and evaluation results
saveRDS(sim_data, 'result/Evaluation/sim-data.rds')
write.csv(sim_embed$mds, 'result/Evaluation/sim-mds-Z.csv', row.names = FALSE)
write.csv(sim_embed$umap_s5, 'result/Evaluation/sim-umap_s-05-Z.csv', row.names = FALSE) # supervised, n_neignbor = 5
write.csv(sim_embed$umap_s10, 'result/Evaluation/sim-umap_s-10-Z.csv', row.names = FALSE) # supervised, n_neignbor = 10
write.csv(sim_embed$umap_s20, 'result/Evaluation/sim-umap_s-20-Z.csv', row.names = FALSE) # supervised, n_neignbor = 20
write.csv(sim_embed$umap_s30, 'result/Evaluation/sim-umap_s-30-Z.csv', row.names = FALSE) # supervised, n_neignbor = 30
write.csv(sim_embed$umap_u5, 'result/Evaluation/sim-umap_u-05-Z.csv', row.names = FALSE) # unsupervised, n_neignbor = 5
write.csv(sim_embed$umap_u10, 'result/Evaluation/sim-umap_u-10-Z.csv', row.names = FALSE) # unsupervised, n_neignbor = 10
write.csv(sim_embed$umap_u20, 'result/Evaluation/sim-umap_u-20-Z.csv', row.names = FALSE) # unsupervised, n_neignbor = 20
write.csv(sim_embed$umap_u30, 'result/Evaluation/sim-umap_u-30-Z.csv', row.names = FALSE) # unsupervised, n_neignbor = 30
write.csv(sim_embed$tsne5, 'result/Evaluation/sim-tsne-05-Z.csv', row.names = FALSE) # perplexity = 5
write.csv(sim_embed$tsne10, 'result/Evaluation/sim-tsne-10-Z.csv', row.names = FALSE) # perplexity = 10
write.csv(sim_embed$tsne20, 'result/Evaluation/sim-tsne-20-Z.csv', row.names = FALSE) # perplexity = 20
write.csv(sim_embed$tsne30, 'result/Evaluation/sim-tsne-30-Z.csv', row.names = FALSE) # perplexity = 30
write.csv(sim_embed$iso5, 'result/Evaluation/sim-iso-05-Z.csv', row.names = FALSE) # no shortest dissimilarities = 5
write.csv(sim_embed$iso10, 'result/Evaluation/sim-iso-10-Z.csv', row.names = FALSE) # no shortest dissimilarities = 10
write.csv(sim_embed$iso20, 'result/Evaluation/sim-iso-20-Z.csv', row.names = FALSE) # no shortest dissimilarities = 20
write.csv(sim_embed$iso30, 'result/Evaluation/sim-iso-30-Z.csv', row.names = FALSE) # no shortest dissimilarities = 30

saveRDS(eval_res2$sim$mds)

# run SMDS using sim_data
library(superMDS)
for(alpha in (1:4)/4){
  res <- TrainSuperMDS(d = sim_data$distmat, y = sim_data$data$Y, alpha = alpha)
  write.csv(res$z, sprintf('result/Evaluation/sim-smds-%.2f-Z.csv', alpha), row.names=FALSE)
  print('step done')
}

# run FMDS using sim_data
library(parallel)
library(pbmcapply)
source('mm.R')
sim_data <- readRDS('result/Evaluation/sim-data.rds')
z0 <- read.csv('result/Evaluation/sim-mds-Z.csv')
res <- pbmclapply(1:4, function(i){
  x <- c(1:4)[i]/5
  return(mm_cmds(nit=50, lambda=x, z0=z0, D=sim_data$distmat, 
                 y = sim_data$data$Y, dataset = 'sim'))
}, mc.cores = 4)


## Algal microbiome
saveRDS(dist1, 'result/Evaluation/alga_1-dist.rds')
saveRDS(dist2, 'result/Evaluation/alga_2-dist.rds')
write.csv(y1, 'result/Evaluation/alga_1-Y.csv', row.names = FALSE)
write.csv(y2, 'result/Evaluation/alga_2-Y.csv', row.names = FALSE)
write.csv(zmds1, 'result/Evaluation/alga_1-mds-Z.csv', row.names = FALSE)
write.csv(zmds2, 'result/Evaluation/alga_2-mds-Z.csv', row.names = FALSE)

write.csv(site2_embed$umap_s5, 'result/Evaluation/alga_2-umap_s-05-Z.csv', row.names = FALSE) # supervised, n_neignbor = 5
write.csv(site2_embed$umap_s10, 'result/Evaluation/alga_2-umap_s-10-Z.csv', row.names = FALSE) # supervised, n_neignbor = 10
write.csv(site2_embed$umap_s20, 'result/Evaluation/alga_2-umap_s-20-Z.csv', row.names = FALSE) # supervised, n_neignbor = 20
write.csv(site2_embed$umap_s30, 'result/Evaluation/alga_2-umap_s-30-Z.csv', row.names = FALSE) # supervised, n_neignbor = 30
write.csv(site2_embed$umap_u5, 'result/Evaluation/alga_2-umap_u-05-Z.csv', row.names = FALSE) # unsupervised, n_neignbor = 5
write.csv(site2_embed$umap_u10, 'result/Evaluation/alga_2-umap_u-10-Z.csv', row.names = FALSE) # unsupervised, n_neignbor = 10
write.csv(site2_embed$umap_u20, 'result/Evaluation/alga_2-umap_u-20-Z.csv', row.names = FALSE) # unsupervised, n_neignbor = 20
write.csv(site2_embed$umap_u30, 'result/Evaluation/alga_2-umap_u-30-Z.csv', row.names = FALSE) # unsupervised, n_neignbor = 30
write.csv(site2_embed$tsne5, 'result/Evaluation/alga_2-tsne-05-Z.csv', row.names = FALSE) # perplexity = 5
write.csv(site2_embed$tsne10, 'result/Evaluation/alga_2-tsne-10-Z.csv', row.names = FALSE) # perplexity = 10
write.csv(site2_embed$tsne7, 'result/Evaluation/alga_2-tsne-07-Z.csv', row.names = FALSE) # perplexity = 7
write.csv(site2_embed$iso5, 'result/Evaluation/alga_2-iso-05-Z.csv', row.names = FALSE) # no shortest dissimilarities = 5
write.csv(site2_embed$iso10, 'result/Evaluation/alga_2-iso-10-Z.csv', row.names = FALSE) # no shortest dissimilarities = 10
write.csv(site2_embed$iso7, 'result/Evaluation/alga_2-iso-07-Z.csv', row.names = FALSE) # no shortest dissimilarities = 7
saveRDS(nn2, 'result/Evaluation/alga_2-nn-data.rds')

# SMDS
library(superMDS)
for(alpha in (1:4)/5){
  res <- TrainSuperMDS(d = distmat1, y = y1, alpha = alpha)
  write.csv(res$z, sprintf('result/Evaluation/alga_1-smds-%.2f-Z.csv', alpha), row.names=FALSE)
  res <- TrainSuperMDS(d = distmat2, y = y2, alpha = alpha)
  write.csv(res$z, sprintf('result/Evaluation/alga_2-smds-%.2f-Z.csv', alpha), row.names=FALSE)
  print('alpha done')
}

# FMDS
source('mm.R')
mm_cmds(nit=50, lambda=0.38, z0=zmds1, D=distmat1, y = y1, dataset = 'alga_1')
mm_cmds(nit=50, lambda=0.38, z0=zmds2, D=distmat2, y = y2, dataset = 'alga_2')
# when lambda exceeds 0.6 it explodes in alga_1 presumably because of the non-negativeness of second order coefficient term.


## cirrhosis
saveRDS(phyl_unifrac_cirrhosis, 'result/Evaluation/cirr-dist.rds')
write.csv(cirr_y, 'result/Evaluation/cirr-Y.csv', row.names = FALSE)
write.csv(cirr_embed$mds, 'result/Evaluation/cirr-mds-Z.csv', row.names = FALSE)

## T2D
saveRDS(phyl_unifrac_cirrhosis, 'result/Evaluation/t2d-dist.rds')
write.csv(cirr_y, 'result/Evaluation/t2d-Y.csv', row.names = FALSE)
write.csv(cirr_embed$mds, 'result/Evaluation/t2d-mds-Z.csv', row.names = FALSE)

## run FMDS for human microbiome
library(parallel)
source('mm.R')
dist <- readRDS('result/Evaluation/t2d-dist.rds')
z0 <- as.matrix(read.csv('result/Evaluation/t2d-mds-Z.csv'))
y <- as.matrix(read.csv('result/Evaluation/t2d-Y.csv'))

res <- mclapply(1:4, function(i){
  x <- c(1:4)[i]/5
  return(mm_cmds(nit=50, lambda=x, z0=z0, D=as.matrix(dist), y = y, dataset = 't2d'))
}, mc.cores = 4)