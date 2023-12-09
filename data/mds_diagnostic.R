# variables dependency on get_dist.R, mm.R, sim_etc.R
source("mm.R")
library(superMDS)

get_stress <- function(dist0, dist){
  res <- sqrt(sum((dist0 - dist)^2) / sum(dist^2))
  return(res)
}

smds_x <- TrainSuperMDS(d = distmat2, y = y2s[,1], alpha = 0.3/(1+0.3))


#### Compare stress
## Simulated superMDS
print(Sys.time())
sim_res_smds <- mclapply(1:4, function(i){
  x <- c(0, 0.3, 0.5, 0.7)[i]
  return(TrainSuperMDS(d = sim_data$distmat, y = sim_data$data$Y, alpha = x/(1+x)))
  print(Sys.time())
  }, mc.cores = 4)
names(sim_res_smds) <- c("lambda0", "lambda0.3", "lambda0.5", "lambda0.7")

print(round(get_stress(sim_data$dist, dist(sim_res_smds$lambda0$z)), 4)) # 0.1309
print(round(get_stress(sim_data$dist, dist(sim_res_smds$lambda0.3$z)), 4)) # 0.3023
print(round(get_stress(sim_data$dist, dist(sim_res_smds$lambda0.5$z)), 4)) # 0.3694
print(round(get_stress(sim_data$dist, dist(sim_res_smds$lambda0.7$z)), 4)) # 0.42


#### Shepard Plot
## simulation superMDS
pdf("result/shepard_sim_smds.pdf", width = 7, height = 8)
par(mfrow = c(2, 2))
plot(sim_data$dist, dist(sim_res_smds$lambda0$z),
     xlab = "original", ylab = "configuration",
     main = "lambda = 0", col = alpha("black", 0.15), pch=16)
plot(sim_data$dist, dist(sim_res_smds$lambda0.3$z),
     xlab = "original", ylab = "configuration",
     main = "lambda = 0.3", col = alpha("black", 0.15), pch=16)
plot(sim_data$dist, dist(sim_res_smds$lambda0.5$z),
     xlab = "original", ylab = "configuration",
     main = "lambda = 0.5", col = alpha("black", 0.15), pch=16)
plot(sim_data$dist, dist(sim_res_smds$lambda0.7$z),
     xlab = "original", ylab = "configuration",
     main = "lambda = 0.7", col = alpha("black", 0.15), pch=16)
dev.off()

## simulation proposed
pdf("result/shepard_sim.pdf", width = 7, height = 8)
par(mfrow = c(2, 2))
plot(sim_data$dist, dist(sim_res$mds),
     xlab = "original", ylab = "configuration",
     main = "lambda = 0", col = alpha("black", 0.15), pch=16)
plot(sim_data$dist, dist(sim_res$proposed$lambda0.3$z),
     xlab = "original", ylab = "configuration",
     main = "lambda = 0.3", col = alpha("black", 0.15), pch=16)
plot(sim_data$dist, dist(sim_res$proposed$lambda0.5$z),
     xlab = "original", ylab = "configuration",
     main = "lambda = 0.5", col = alpha("black", 0.15), pch=16)
plot(sim_data$dist, dist(sim_res$proposed$lambda0.3$z),
     xlab = "original", ylab = "configuration",
     main = "lambda = 0.7", col = alpha("black", 0.15), pch=16)
dev.off()

#### Spearman correlation from Shepard
## simulation superMDS
print(round(
  cov(sim_data$dist, dist(sim_res_smds$lambda0$z)) / 
    (sd(sim_data$dist) * sd(dist(sim_res_smds$lambda0$z))), 
  4)) # 0.9649
print(round(
  cov(sim_data$dist, dist(sim_res_smds$lambda0.3$z)) / 
    (sd(sim_data$dist) * sd(dist(sim_res_smds$lambda0.3$z))), 
  4)) # 0.8139
print(round(
  cov(sim_data$dist, dist(sim_res_smds$lambda0.5$z)) / 
    (sd(sim_data$dist) * sd(dist(sim_res_smds$lambda0.5$z))), 
  4)) # 0.7479
print(round(
  cov(sim_data$dist, dist(sim_res_smds$lambda0.7$z)) / 
    (sd(sim_data$dist) * sd(dist(sim_res_smds$lambda0.7$z))), 
  4)) # 0.6981

## simulation  proposed MDS
print(round(
  cov(sim_data$dist, dist(sim_res$mds)) / 
    (sd(sim_data$dist) * sd(dist(sim_res$mds))), 
  4)) # 0.9594
print(round(
  cov(sim_data$dist, dist(sim_res$proposed$lambda0.3$z)) / 
    (sd(sim_data$dist) * sd(dist(sim_res$proposed$lambda0.3$z))), 
  4)) # 0.9567
print(round(
  cov(sim_data$dist, dist(sim_res$proposed$lambda0.5$z)) / 
    (sd(sim_data$dist) * sd(dist(sim_res$proposed$lambda0.5$z))), 
  4)) # 0.9317
print(round(
  cov(sim_data$dist, dist(sim_res$proposed$lambda0.7$z)) / 
    (sd(sim_data$dist) * sd(dist(sim_res$proposed$lambda0.7$z))), 
  4)) # 0.9013


# Some definitions switches x and y axis
pdf("shepard.pdf", width = 10, height = 3)
par(mfrow = c(1,3))
plot(dist2, dist(zmds2),
     xlab = "original distance", ylab = "configuration distance",
     main = "Pure MDS")
mtext(paste(
  "Stress1 = ", 
  round(sqrt(sum((dist2 - dist(zmds2))^2) /
               sum((dist(zmds2))^2)), 4)), side=3)
plot(dist2, dist(obmm_x$z),
     xlab = "original distance", ylab = "configuration distance",
     main = expression(paste(lambda, "= 0.3")))
mtext(paste(
  "Stress1 = ", 
  round(sqrt(sum((dist2 - dist(obmm_x$z))^2) /
               sum((dist(obmm_x$z))^2)), 4)), side=3)
plot(dist2, dist(smds_x$z),
     xlab = "original distance", ylab = "configuration distance",
     main = expression(paste("SMDS ", alpha, "= 0.231")))
mtext(paste(
  "Stress1 = ", 
  round(sqrt(sum((dist2 - dist(smds_x$z))^2) /
               sum((dist(smds_x$z))^2)), 4)), side=3)
dev.off()
## Stress-1: Scaled Stress
## For its definition, see, for example,
## Andreas Buja, Deborah F Swayne, Michael L Littman, Nathaniel Dean, Heike Hofmann & Lisha Chen (2008) 
## Data Visualization With Multidimensional Scaling, 
## Journal of Computational and Graphical Statistics, 17:2, 444-472
## (it is in our shared google drive)

### Stress-Per-Point (SPP)
## This is usually used for outlier detection 
## (detects ones with extraordinarily large Stress)
spp <- function(D = distmat2, z = obmm_x$z){
  N <- nrow(D)
  spp <- apply((D - get_dist_mat(z))^2, 1, function(x){sum(x)/(N-1)})
  return(spp)
}
pdf("spp_bubble.pdf", width = 5, height = 4)
ggplot(data = data.frame(obmm_x$z, y = as.factor(y2s[,1]), spp = spp()),
       aes(x = Axis.1, y = Axis.2, size = spp,
           color = y)) +
  geom_point() +
  labs(title="Stress-Per-Point Bubble Plot")
dev.off()
