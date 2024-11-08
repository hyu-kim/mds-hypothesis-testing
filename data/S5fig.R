source('fig_util.R')
library(scales)
# custom plot option
myplot <- function(x, y, t, p, ...){
  plot(x, y, xlab = "", ylab = "", xaxt = 'n', yaxt ='n',
       main=sprintf("%s, %s",t, p),
       col = alpha("black", 0.3),
       pch=16, cex=0.4, cex.axis=0.8, cex.main=1,
       xlim = c(0, 10), ylim = c(0, 10),
       tcl=-0.2, lwd=0.5, ...
  )
  axis(side=1, lwd=0.5, tck=-0.05)
  axis(side=2, lwd=0.5, tck=-0.05, las=2)
}


plot_fcor_set <- function(dataset='alga_2', method='fmds', method_f='FMDS', params_v=c(0.5,1)){
  x_dist <- readRDS(sprintf('result/Evaluation/%s-dist.rds', dataset))
  y_orig <- as.matrix(read.csv(sprintf('result/Evaluation/%s-Y.csv', dataset)))
  len <- length(params_v)
  
  for(p in params_v){
    if(p<=1) {
      z_embed <- as.matrix(read.csv(sprintf('result/Evaluation/%s-%s-%.2f-Z.csv',
                                            dataset,method,p)))
      }
    else {
      z_embed <- as.matrix(read.csv(sprintf('result/Evaluation/%s-%s-%02d-Z.csv',
                                            dataset,method,p)))
      }
      
    res <- rp_eval2(dm=as.matrix(x_dist), y_orig=y_orig, z_emb=z_embed)
    myplot(res$F0_perm, res$Femb_perm, method_f, p)
    points(res$F0, res$Femb, pch=4, col='red')
    print(paste(res$F0, res$Femb))
  }
  
  while(len<4){
    plot.new()
    len <- len+1
  }
}

method_v <- c('mds','nn','fmds','smds','umap_s','umap_u','tsne','iso')
method_formal_v <- c('MDS','NN','FMDS','SMDS','UMAP-S','UMAP-U','t-SNE','Isomap')
params_v1 <- (1:4)/5 # fmds, smds
params_v2 <- c(5,10,20,30) # umap
params_v3 <- c(5,7,10) # tsne, isomap


pdf("result/Fig_S5.pdf", width = 5, height = 8)

par(mfrow = c(7,4), mar = c(1.5,2.5,0.8,0.1), mgp = c(0,0.3,0), lwd=0.75)

for(i in seq(8)){
  m <- method_v[i]
  m_f <- method_formal_v[i]
  x_dist <- readRDS('result/Evaluation/alga_2-dist.rds')
  y_orig <- as.matrix(read.csv('result/Evaluation/alga_2-Y.csv'))
  
  if(m=='mds'){
    z_embed <- as.matrix(read.csv('result/Evaluation/alga_2-mds-Z.csv'))
    res <- rp_eval2(dm=as.matrix(x_dist), y_orig=y_orig, z_emb=z_embed)
    myplot(res$F0_perm, res$Femb_perm, 'MDS', '')
    points(res$F0, res$Femb, pch=4, col='red')
    plot.new()
    print(m)
    next
  }
  else if(m=='nn'){
    alga_2_nn <- readRDS('result/Evaluation/alga_2-nn-data.rds')
    z_embed <- alga_2_nn[,2:33]
    res <- rp_eval2(dm=as.matrix(x_dist), y_orig=y_orig, z_emb=z_embed)
    myplot(res$F0_perm, res$Femb_perm, 'NN', '')
    points(res$F0, res$Femb, pch=4, col='red')
    plot.new()
    print(m)
    next
  }
  
  if(m=='fmds' | m=='smds') params_v <- params_v1
  else if(m=='umap_s' | m=='umap_u') params_v <- params_v2
  else params_v <- params_v3
  
  plot_fcor_set(method=m, method_f=m_f, params_v=params_v)
  print(m)
}

dev.off()