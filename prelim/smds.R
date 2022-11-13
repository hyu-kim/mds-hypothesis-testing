# install.packages("superMDS")
source("pcoa+permanova.R")
library(superMDS)

#### Get dissimilarity matrix
site1 = subset_samples(ps, Site=="1") # subset
site2 = subset_samples(ps, Site=="2") # subset

dist1 = phyloseq::distance(site1, method = "unifrac")
dist2 = phyloseq::distance(site2, method = "unifrac")

distmat1 = distmat2 = matrix(0, nrow = 36, ncol = 36)
pnt = 0
for(i in 1:35){
  len = 36-i
  distmat1[(i+1):36, i] = dist1[(pnt+1):(pnt+len)]
  distmat2[(i+1):36, i] = dist2[(pnt+1):(pnt+len)]
  pnt = pnt+len
}
distmat1 = distmat1 + t(distmat1)
distmat2 = distmat2 + t(distmat2)

y1 <- ifelse(site1@sam_data$Treatment == "Pt +", 1, 2)
y2 <- ifelse(site2@sam_data$Treatment == "Pt +", 1, 2)

#### SMDS to our data
smds_site1 <- list(
  a0 = TrainSuperMDS(d = distmat1, y = y1, alpha = 0.0),
  a5 = TrainSuperMDS(d = distmat1, y = y1, alpha = 0.5),
  a1 = TrainSuperMDS(d = distmat1, y = y1, alpha = 1)
)
smds_site2 <- list(
  a0 = TrainSuperMDS(d = distmat2, y = y2, alpha = 0.0),
  a5 = TrainSuperMDS(d = distmat2, y = y2, alpha = 0.5),
  a1 = TrainSuperMDS(d = distmat2, y = y2, alpha = 1)
)

# pdf("SMDS_res.pdf", width = 7.5, height = 5)
par(mfrow = c(2, 3), mar = c(2,2,2,1))
plot(smds_site1$a0$z, col = ifelse(y1 == 1, 1, 2), 
     pch = ifelse(y1 == 1, 16, 17), cex = 1.5,
     xlab = "", ylab = "", main = "Site 1, alpha = 0")
plot(smds_site1$a5$z, col = ifelse(y1 == 1, 1, 2), 
     pch = ifelse(y1 == 1, 16, 17), cex = 1.5,
     xlab = "", ylab = "", main = "Site 1, alpha = 0.5")
plot(smds_site1$a1$z, col = ifelse(y1 == 1, 1, 2), 
     pch = ifelse(y1 == 1, 16, 17), cex = 1.5,
     xlab = "", ylab = "", main = "Site 1, alpha = 1")
plot(smds_site2$a0$z, col = ifelse(y1 == 1, 1, 2), 
     pch = ifelse(y1 == 1, 16, 17), cex = 1.5,
     xlab = "", ylab = "", main = "Site 2, alpha = 0")
plot(smds_site2$a5$z, col = ifelse(y1 == 1, 1, 2), 
     pch = ifelse(y1 == 1, 16, 17), cex = 1.5,
     xlab = "", ylab = "", main = "Site 2, alpha = 0.5")
plot(smds_site2$a1$z, col = ifelse(y1 == 1, 1, 2), 
     pch = ifelse(y1 == 1, 16, 17), cex = 1.5,
     xlab = "", ylab = "", main = "Site 2, alpha = 1")
legend("bottomright", c("Pt +", "Pt -"),
       col = c(1, 2), pch = c(16, 17))
# dev.off()

#### Simulation
set.seed(10101)
library(MASS)
x <- mvrnorm(100, rep(0, 10), diag(1, 10, 10))
y_lab <- c(rep(1, 50), rep(2, 50))
x_dist <- dist(x)
x_distmat <- matrix(0, nrow = 100, ncol = 100)
pnt = 0
for(i in 1:99){
  len = 100-i
  x_distmat[(i+1):100, i] = x_dist[(pnt+1):(pnt+len)]
  pnt = pnt+len
}
x_distmat = x_distmat + t(x_distmat)
x_smds <- list(
  a0 = TrainSuperMDS(d = x_distmat, y = y_lab, alpha = 0.0),
  a2 = TrainSuperMDS(d = x_distmat, y = y_lab, alpha = 0.2),
  a5 = TrainSuperMDS(d = x_distmat, y = y_lab, alpha = 0.5),
  a1 = TrainSuperMDS(d = x_distmat, y = y_lab, alpha = 1)
)

pdf("sim_smds.pdf", height = 5, width = 5)
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(x_smds$a0$z, col = y_lab, 
     pch = ifelse(y_lab == 1, 16, 17), cex = 1.5,
     xlab = "", ylab = "", main = "alpha = 0")
plot(x_smds$a2$z, col = y_lab, 
     pch = ifelse(y_lab == 1, 16, 17), cex = 1.5,
     xlab = "", ylab = "", main = "alpha = 0.2")
plot(x_smds$a5$z, col = y_lab, 
     pch = ifelse(y_lab == 1, 16, 17), cex = 1.5,
     xlab = "", ylab = "", main = "alpha = 0.5")
plot(x_smds$a1$z, col = y_lab, 
     pch = ifelse(y_lab == 1, 16, 17), cex = 1.5,
     xlab = "", ylab = "", main = "alpha = 1")
dev.off()
