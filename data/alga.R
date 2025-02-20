library(MASS)
library(parallel)
library(pbmcapply)
source('mm.R')

distmat2 <- as.matrix(readRDS('result/Evaluation/alga_2-dist.rds'))
y2 <- as.matrix(read.csv('result/Evaluation/alga_2-Y.csv'))
zmds2 <- read.csv('result/Evaluation/alga_2-mds-Z.csv')

print('importing done. Begins iteration..')
res <- pbmclapply(1:8, function(i){
  x <- c(13:20)[i]/20
  return(mm_cmds(nit=50, 
                 lambda=x, 
                 z0=zmds2, 
                 D=as.matrix(distmat2), 
                 y = y2,
                 dataset = 'alga_2')
  )
}, mc.cores = 4)

res <- mm_cmds(nit=200, lambda=0.25, z0=zmds2, D=as.matrix(dist2), y = y2, dataset = 'alga_2')

# SuperMDS
library(superMDS)
for(alpha in (1:4)/5){
  res <- TrainSuperMDS(d = distmat2, y = y2, alpha = alpha)
  write.csv(res$z, sprintf('result/Evaluation/alga_2-smds-%.2f-Z.csv', alpha), row.names=FALSE)
  print('alpha done')
}

# Scheduled lambda method
model_loess <- readRDS('result/HyperparameterStudy/Nonparametric/model_loess.Rds')
res <- mm_cmds2(nit=200, z0=zmds2, D=as.matrix(distmat2), y = y2, dataset = 'alga_2', model_lambda=model_loess)