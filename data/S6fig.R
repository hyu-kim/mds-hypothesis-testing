pdf("figures/Fig_S6.pdf", width = 6.5, height = 6.5)

par(mfrow = c(3, 3), mar = c(3,3,1,1), mgp = c(2,1,0))

for(r in 1:3){
  data_app <- read.csv('result/Multiclass/sim4d-data.csv')
  y_app <- as.numeric(as.matrix(read.csv('result/Multiclass/sim4d-Y.csv')))
  
  plot(data_app$X1, data_app$X4,
       xlab = "X1", ylab = "X4", 
       pch = c(16,17,15)[y_app], col = c('red', 'blue', 'green4')[y_app])
  plot(data_app$X2, data_app$X4,
       xlab = "X2", ylab = "X4", 
       pch = c(16,17,15)[y_app], col = c('red', 'blue', 'green4')[y_app])
  plot(data_app$X3, data_app$X4,
       xlab = "X3", ylab = "X4", 
       pch = c(16,17,15)[y_app], col = c('red', 'blue', 'green4')[y_app])
  
  plot(data_app$X1, data_app$X3,
       xlab = "X1", ylab = "X3",
       pch = c(16,17,15)[y_app], col = c('red', 'blue', 'green4')[y_app])
  plot(data_app$X2, data_app$X3,
       xlab = "X2", ylab = "X3",
       pch = c(16,17,15)[y_app], col = c('red', 'blue', 'green4')[y_app])
  plot.new()
  
  plot(data_app$X1, data_app$X2,
       xlab = "X1", ylab = "X2",
       pch = c(16,17,15)[y_app], col = c('red', 'blue', 'green4')[y_app])
  plot.new()
  plot.new()
}

dev.off()