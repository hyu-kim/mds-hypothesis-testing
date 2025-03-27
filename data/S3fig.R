source('fig_util.R')
library(scales)

# custom plot option
myplot <- function(x, y, t, p, ...){
  plot(x, y, xlab = "", ylab = "", xaxt = 'n', yaxt ='n',
       main=sprintf("%s, %s",t, p),
       col = alpha("black", 0.3),
       pch=16, cex=0.6, cex.axis=0.8, cex.main=0.8,
       xlim = c(0, 9), ylim = c(0, 9),
       tcl=-0.2, lwd=0.5, ...
  )
  axis(side=1, lwd=0.5, tck=-0.05)
  axis(side=2, lwd=0.5, tck=-0.05, las=2)
}


plot_fcor_set <- function(dataset='sim_1', method='fmds', method_f='FMDS', params_v=c(0.5,1)){
  data_df <- as.matrix(read.csv(sprintf('Data/Simulated/%s-data.csv', dataset)))
  y_orig <- as.matrix(read.csv(sprintf('Data/Simulated/%s-Y.csv', dataset)))
  
  for(p in params_v){
    if(p<=1) {
      z_embed <- as.matrix(read.csv(sprintf('Data/Simulated/%s/%s-%s-%.2f-Z.csv',
                                            method_f,dataset,method,p)))
      }
    else {
      z_embed <- as.matrix(read.csv(sprintf('Data/Simulated/%s/%s-%s-%02d-Z.csv',
                                            method_f,dataset,method,p)))
      }
      
    res <- rp_eval2(dm=as.matrix(dist(data_df)), y_orig=y_orig, z_emb=z_embed)
    myplot(res$F0_perm, res$Femb_perm, method_f, p)
    points(res$F0, res$Femb, pch=4, col='red')
    print(paste(res$F0, res$Femb))
  }
}

method_v <- c('mds','fmds','smds','umap_s','umap_u','tsne','iso')
method_formal_v <- c('MDS','F-MDS','superMDS','UMAP-S','UMAP-U','t-SNE','Isomap')
params_v1 <- c(1:4)/5
params_v2 <- c(5,10,20,30)


pdf("figures/Fig_S3.pdf", width = 5, height = 8)

par(mfrow = c(7,4), mar = c(1.5,2.5,0.8,0.1), mgp = c(0,0.3,0), lwd=0.75)

for(i in seq(7)){
  m <- method_v[i]
  m_f <- method_formal_v[i]
  
  if(m=='mds'){
    data_df <- as.matrix(read.csv('result/Evaluation/sim_1-data.csv'))
    z_embed <- as.matrix(read.csv('result/Evaluation/sim_1-mds-Z.csv'))
    y_orig <- as.matrix(read.csv('result/Evaluation/sim_1-Y.csv'))
    res <- rp_eval2(dm=as.matrix(dist(data_df)), y_orig=y_orig, z_emb=z_embed)
    myplot(res$F0_perm, res$Femb_perm, 'MDS', '')
    points(res$F0, res$Femb, pch=4, col='red')
    plot.new()
    plot.new()
    plot.new()
    print(m)
    next
  }
  
  if(m=='fmds' | m=='smds') params_v <- params_v1
  else params_v <- params_v2
  
  plot_fcor_set(method=m, method_f=m_f, params_v=params_v)
  print(m)
}

dev.off()