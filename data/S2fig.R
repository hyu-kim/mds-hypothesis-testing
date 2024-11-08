# custom plot option
myplot <- function(x, y, t, p, ...){
  plot(x, y, xlab = "", ylab = "", xaxt = 'n', yaxt ='n',
       main=sprintf("%s, %s",t, p),
       col = alpha("black", 0.15), 
       pch=16, cex=0.3, cex.axis=0.8, cex.main=0.8,
       tcl=-0.2, lwd=0.5, ...
  )
  axis(side=1, lwd=0.5, tck=-0.05)
  axis(side=2, lwd=0.5, tck=-0.05, las=2)
}

plot_shepard_set <- function(dataset='sim_1', method='fmds', method_f='FMDS', params_v=c(0.5,1)){
  x_dist <- get_dist_mat(read.csv(sprintf('result/Evaluation/%s-data.csv', dataset)))
  
  for(p in params_v){
    if(p<=1) {
      z_dist <- get_dist_mat(
        read.csv(sprintf('result/Evaluation/%s-%s-%.2f-Z.csv',dataset,method,p)))
    }
    else {z_dist <- get_dist_mat(
      read.csv(sprintf('result/Evaluation/%s-%s-%02d-Z.csv',dataset,method,p)))}
    
    myplot(x_dist, z_dist, method_f, p)
  }
}


method_v <- c('mds','fmds','smds','umap_s','umap_u','tsne','iso')
method_formal_v <- c('MDS','FMDS','SMDS','UMAP-S','UMAP-U','t-SNE','Isomap')
params_v1 <- c(1:4)/5
params_v2 <- c(5,10,20,30)


pdf("figures/Fig_S2.pdf", width = 5, height = 8)

# par(mfrow = c(7,4), mar = c(3,3,1,1), mgp = c(2,1,0))
par(mfrow = c(7,4), mar = c(1.5,2.5,0.8,0.1), mgp = c(0,0.3,0), lwd=0.75)

for(i in seq(7)){
  m <- method_v[i]
  m_f <- method_formal_v[i]
  if(m=='mds'){
    x_dist <- get_dist_mat(read.csv('result/Evaluation/sim_1-data.csv'))
    z_dist <- get_dist_mat(read.csv('result/Evaluation/sim_1-mds-Z.csv'))
    myplot(x_dist, z_dist, 'MDS', '')
    plot.new()
    plot.new()
    plot.new()
    print(m)
    next
  }
  
  if(m=='fmds' | m=='smds') params_v <- params_v1
  else params_v <- params_v2
  
  plot_shepard_set(method=m, method_f=m_f, params_v=params_v)
  print(m)
}

dev.off()