library(ape)

features = read.csv('result/nn/features_out.csv', header=FALSE)
labels = read.csv('result/nn/labels_out.csv', header=FALSE)

simclr_dist <- dist(features)

features_pcoa <- pcoa(simclr_dist)

postscript('result/nn/pcoa_nn.eps', width = 6, height = 4)
plot(features_pcoa$vectors[,1:2], col=factor(lab))
dev.off()

dist_all = phyloseq::distance(ps, method = "unifrac", weighted = T)

postscript('result/nn/shepard_nn.eps', width = 4, height = 4)
plot(dist_all/max(dist_all), simclr_dist/max(simclr_dist), 
     xlab = "original distance", ylab = "features distance",
     pch='.')
mtext(paste(
  "Stress1 = ", 
  round(sqrt(sum((dist_all/max(dist_all) - simclr_dist/max(simclr_dist))^2) /
               sum((simclr_dist/max(simclr_dist))^2)), 4)), 
  side=3)

mtext(paste(
  "Correlation = ", 
  round(cor(dist_all, simclr_dist), 
        4)
  ), 
  side=3, line=1)
dev.off()