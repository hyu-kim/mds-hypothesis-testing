#####TTD
#0. Do with real data labels
#### Site 1
res1 <- 
  list(
    lambda0.0 = mm_cmds(nit=15, lambda=0.0, z0=zmds1, D=distmat1, y=y1),
    lambda0.1 = mm_cmds(nit=15, lambda=0.1, z0=zmds1, D=distmat1, y=y1),
    lambda0.3 = mm_cmds(nit=15, lambda=0.3, z0=zmds1, D=distmat1, y=y1),
    lambda0.5 = mm_cmds(nit=15, lambda=0.5, z0=zmds1, D=distmat1, y=y1),
    lambda10= mm_cmds(nit=15, lambda=010, z0=zmds1, D=distmat1, y=y1)
  )

## Desired result (PERMANOVA)
get_p(trt = y1, d = distmat1)$ratio # pseudo-F = 7.402205
get_p(trt = y1, d = distmat1)$p # p = 0

## Pure MDS
get_p(trt = y1, mat = zmds1)$ratio # pseudo-F = 11.984
get_p(trt = y1, mat = zmds1)$p # p = 0

## Proposed MDS
res1$lambda0.3$F_z # pseudo-F = 12.36596
get_p(trt = y1, mat = res1$lambda0.3$z)$p # p = 0

## Configurations
pdf("data/result/config_site1.pdf",
    width = 9, height = 5)
par(mfrow = c(1,2))
plot(zmds1, main = "Pure MDS", col = y1)
text(x = (min(zmds1[,1])+max(zmds1[,1]))/2, y = max(zmds1[,2]),
     paste("F_z = ", round(get_p(trt = y1, mat = zmds1)$ratio, 3)))
text(x = (min(zmds1[,1])+max(zmds1[,1]))/2, y = max(zmds1[,2]) - 0.008,
     paste("p = ", get_p(trt = y1, mat = zmds1)$p))
mtext(paste(
  "Stress1 = ", 
  round(sqrt(sum((dist1 - dist(zmds1))^2) /
               sum((dist(zmds1))^2)), 4)), side=3)

plot(res1$lambda0.3$z, main = "lambda = 0.3", col = y1)
text(x = (min(res1$lambda0.3$z[,1])+max(res1$lambda0.3$z[,1]))/2, 
     y = max(res1$lambda0.3$z[,2]),
     paste("F_z = ", round(res1$lambda0.3$F_z, 3)))
text(x = (min(res1$lambda0.3$z[,1])+max(res1$lambda0.3$z[,1]))/2, 
     y = max(res1$lambda0.3$z[,2]) - 0.008,
     paste("p = ", get_p(trt = y1, mat = res1$lambda0.3$z)$p))
mtext(paste(
  "Stress1 = ", 
  round(sqrt(sum((dist1 - dist(res1$lambda0.3$z))^2) /
               sum((dist(res1$lambda0.3$z))^2)), 4)), side=3)
dev.off()

## Shepard
pdf("data/result/shepard_site1.pdf",
    width = 9, height = 5)
par(mfrow = c(1, 2))
plot(dist1, dist(zmds1),
     xlab = "original distance", ylab = "configuration distance",
     main = "Pure MDS")
mtext(paste(
  "Stress1 = ", 
  round(sqrt(sum((dist1 - dist(zmds1))^2) /
               sum((dist(zmds1))^2)), 4)), side=3)
plot(dist1, dist(res1$lambda0.3$z),
     xlab = "original distance", ylab = "configuration distance",
     main = expression(paste(lambda, "= 0.3")))
mtext(paste(
  "Stress1 = ", 
  round(sqrt(sum((dist1 - dist(res1$lambda0.3$z))^2) /
               sum((dist(res1$lambda0.3$z))^2)), 4)), side=3)
dev.off()


#### Site 2
res2 <- 
  list(
    lambda0.0 = mm_cmds(nit=15, lambda=0.0, z0=zmds2, D=distmat2, y=y2),
    lambda0.1 = mm_cmds(nit=15, lambda=0.1, z0=zmds2, D=distmat2, y=y2),
    lambda0.3 = mm_cmds(nit=15, lambda=0.3, z0=zmds2, D=distmat2, y=y2),
    lambda0.5 = mm_cmds(nit=15, lambda=0.5, z0=zmds2, D=distmat2, y=y2),
    lambda10 = mm_cmds(nit=15, lambda=10, z0=zmds2, D=distmat2, y=y2)
  )

## Desired result (PERMANOVA)
get_p(trt = y2, d = distmat2)$ratio # pseudo-F = 1.925012
get_p(trt = y2, d = distmat2)$p # p = 0.076

## Pure MDS
get_p(trt = y2, mat = zmds2)$ratio # pseudo-F = 0.733862
get_p(trt = y2, mat = zmds2)$p # p = 0.47

## Proposed MDS
res2$lambda0.3$F_z # pseudo-F = 2.403895
get_p(trt = y2, mat = res2$lambda0.3$z)$p # p = 0.108

## Configurations
pdf("data/result/config_site2.pdf",
    width = 9, height = 5)
par(mfrow = c(1,2))
plot(zmds2, main = "Pure MDS", col = y2)
text(x = (min(zmds2[,1])+max(zmds2[,1]))/2, y = max(zmds2[,2]),
     paste("F_z = ", round(get_p(trt = y2, mat = zmds2)$ratio, 3)))
text(x = (min(zmds2[,1])+max(zmds2[,1]))/2, y = max(zmds2[,2]) - 0.008,
     paste("p = ", get_p(trt = y2, mat = zmds2)$p))
mtext(paste(
  "Stress1 = ", 
  round(sqrt(sum((dist2 - dist(zmds2))^2) /
               sum((dist(zmds2))^2)), 4)), side=3)

plot(res2$lambda0.3$z, main = "lambda = 0.3", col = y2)
text(x = (min(res2$lambda0.3$z[,1])+max(res2$lambda0.3$z[,1]))/2, 
     y = max(res2$lambda0.3$z[,2]),
     paste("F_z = ", round(res2$lambda0.3$F_z, 3)))
text(x = (min(res2$lambda0.3$z[,1])+max(res2$lambda0.3$z[,1]))/2, 
     y = max(res2$lambda0.3$z[,2]) - 0.008,
     paste("p = ", get_p(trt = y2, mat = res2$lambda0.3$z)$p))
mtext(paste(
  "Stress1 = ", 
  round(sqrt(sum((dist2 - dist(res2$lambda0.3$z))^2) /
               sum((dist(res2$lambda0.3$z))^2)), 4)), side=3)
dev.off()

## Shepard
pdf("data/result/shepard_site2.pdf",
    width = 9, height = 5)
par(mfrow = c(1, 2))
plot(dist2, dist(zmds2),
     xlab = "original distance", ylab = "configuration distance",
     main = "Pure MDS")
mtext(paste(
  "Stress1 = ", 
  round(sqrt(sum((dist2 - dist(zmds2))^2) /
               sum((dist(zmds2))^2)), 4)), side=3)
plot(dist2, dist(res2$lambda0.3$z),
     xlab = "original distance", ylab = "configuration distance",
     main = expression(paste(lambda, "= 0.3")))
mtext(paste(
  "Stress1 = ", 
  round(sqrt(sum((dist2 - dist(res2$lambda0.3$z))^2) /
               sum((dist(res2$lambda0.3$z))^2)), 4)), side=3)
dev.off()

#1. Check no symm-matrix issues (summing twice, etc)
#2. Try with permutation on original data (previous to unifrac)

#3. Simulation
library(MASS)
library(plotly)
set.seed(100)
sim_data <- list(
  data = rbind(
  mvrnorm(50, c(0,0,0), Sigma = matrix(c(3,0,0,0,3,0,0,0,1), nrow=3)),
  mvrnorm(50, c(0,0,1), Sigma = matrix(c(3,0,0,0,3,0,0,0,1), nrow=3)))
)
sim_data$data <- data.frame(sim_data$data)
sim_data$data$Y <- c(rep(1, 50), rep(2, 50))
par(mfrow = c(1,1))
plot_ly(sim_data$data, x = ~X1, y = ~X2, z = ~X3, 
        color = ~as.factor(Y), colors = c("red", "black"))
sim_data$distmat <- as.matrix(dist(sim_data$data[, 1:3]))
sim_data$dist <- dist(sim_data$data[, 1:3])

pdf("data/result/sim_data_2d.pdf",
    width = 9, height = 7)
par(mfrow = c(2, 2), mar = c(4,4,2,1))
plot(sim_data$data$X1, sim_data$data$X2, col = sim_data$data$Y,
     xlab = "X1", ylab = "X2")
plot.new()
plot(sim_data$data$X1, sim_data$data$X3, col = sim_data$data$Y,
     xlab = "X1", ylab = "X3")
plot(sim_data$data$X2, sim_data$data$X3, col = sim_data$data$Y,
     xlab = "X2", ylab = "X3")
dev.off()



## Desired result (PERMANOVA)
set.seed(1000)
sim_res$permanova <- get_p(trt = sim_data$data$Y, d = sim_data$distmat)
sim_res$permanova$ratio # pseudo-F = 5.402271
sim_res$permanova$p # p = 0.005

### Pure MDS
sim_res <- list(
  mds = cmdscale(dist(sim_data$data[,1:3]), k = 2))
sim_res$mds <- data.frame(sim_res$mds)
plot(sim_res$mds, col = sim_data$data$Y)
set.seed(100)
sim_res$mds_perm <- get_p(trt = sim_data$data$Y, mat = sim_res$mds)
sim_res$mds_perm$ratio # pseudo-F = 0.07825293
sim_res$mds_perm$p # p = 0.914


## Proposed MDS
library(parallel)
detectCores()
RNGkind("L'Ecuyer-CMRG")
set.seed(10000)
# sim_res$proposed <- 
#   lapply(c(0.1, 0.3, 0.5, 10),
#          FUN = function(x){
#            mm_cmds(nit=15, lambda=x, z0=sim_res$mds, D=sim_data$distmat, 
#                    y = sim_data$data$Y)
#            })
print(Sys.time())
sim_res$proposed <- mclapply(1:3, function(i){
  x <- c(0.3, 0.5, 0.7)[i]
  return(mm_cmds(nit=10, lambda=x, 
                 z0=sim_res$mds, D=sim_data$distmat, 
                 y = sim_data$data$Y))
  print(Sys.time())
  }, mc.cores = 4)
names(sim_res$proposed) <- c("lambda0.3", "lambda0.5", "lambda0.7")

(sim_res$proposed_perm <- c(
  sim_res$proposed$lambda0.3$F_z, # pseudo-F = 2.114592
  get_p(trt = sim_data$data$Y, mat = sim_res$proposed$lambda0.3$z)$p, # p = 0.118
  
  sim_res$proposed$lambda0.5$F_z, # pseudo-F = 5.937904
  get_p(trt = sim_data$data$Y, mat = sim_res$proposed$lambda0.5$z)$p, # p = 0.003
  
  sim_res$proposed$lambda0.7$F_z, # pseudo-F = 5.849832
  get_p(trt = sim_data$data$Y, mat = sim_res$proposed$lambda0.7$z)$p # p = 0.006
  # p values, especially, are subject to change
))



## Configurations
pdf("data/result/config_sim.pdf",
    width = 9, height = 9)
par(mfrow = c(2,2), mar = c(2,2,3,1))
plot(sim_res$mds, main = "Pure MDS", col = sim_data$data$Y)
text(x = (min(sim_res$mds[,1])+
            max(sim_res$mds[,1]))/2 -3, 
     y = max(sim_res$mds[,2]),
     paste("F_z = ", round(sim_res$mds_perm$ratio, 3)))
text(x = (min(sim_res$mds[,1])+
            max(sim_res$mds[,1]))/2 -3, 
     y = max(sim_res$mds[,2]) - 0.5,
     paste("p = ", sim_res$mds_perm$p))
mtext(paste(
  "Stress1 = ", 
  round(sqrt(sum((sim_data$dist - dist(sim_res$mds))^2) /
               sum((dist(sim_res$mds))^2)), 4)), side=3)

plot(sim_res$proposed$lambda0.3$z, main = "lambda = 0.3", col = sim_data$data$Y)
text(x = (min(sim_res$proposed$lambda0.3$z[,1])+
            max(sim_res$proposed$lambda0.3$z[,1]))/2 -3, 
     y = max(sim_res$proposed$lambda0.3$z[,2]),
     paste("F_z = ", round(sim_res$proposed_perm[1], 3)))
text(x = (min(sim_res$proposed$lambda0.3$z[,1])+
            max(sim_res$proposed$lambda0.3$z[,1]))/2 -3, 
     y = max(sim_res$proposed$lambda0.3$z[,2]) - 0.5,
     paste("p = ", sim_res$proposed_perm[2]))
mtext(paste(
  "Stress1 = ", 
  round(sqrt(sum((sim_data$dist - dist(sim_res$proposed$lambda0.3$z))^2) /
               sum((dist(sim_res$proposed$lambda0.3$z))^2)), 4)), side=3)

plot(sim_res$proposed$lambda0.5$z, main = "lambda = 0.5", col = sim_data$data$Y)
text(x = (min(sim_res$proposed$lambda0.5$z[,1])+
            max(sim_res$proposed$lambda0.5$z[,1]))/2 -3, 
     y = max(sim_res$proposed$lambda0.5$z[,2]),
     paste("F_z = ", round(sim_res$proposed_perm[3], 3)))
text(x = (min(sim_res$proposed$lambda0.5$z[,1])+
            max(sim_res$proposed$lambda0.5$z[,1]))/2 -3, 
     y = max(sim_res$proposed$lambda0.5$z[,2]) - 0.5,
     paste("p = ", sim_res$proposed_perm[4]))
mtext(paste(
  "Stress1 = ", 
  round(sqrt(sum((sim_data$dist - dist(sim_res$proposed$lambda0.5$z))^2) /
               sum((dist(sim_res$proposed$lambda0.5$z))^2)), 4)), side=3)

plot(sim_res$proposed$lambda0.7$z, main = "lambda = 0.7", col = sim_data$data$Y)
text(x = (min(sim_res$proposed$lambda0.7$z[,1])+
            max(sim_res$proposed$lambda0.7$z[,1]))/2 -3, 
     y = max(sim_res$proposed$lambda0.7$z[,2]),
     paste("F_z = ", round(sim_res$proposed_perm[5], 3)))
text(x = (min(sim_res$proposed$lambda0.7$z[,1])+
            max(sim_res$proposed$lambda0.7$z[,1]))/2 -3, 
     y = max(sim_res$proposed$lambda0.7$z[,2]) - 0.5,
     paste("p = ", sim_res$proposed_perm[6]))
mtext(paste(
  "Stress1 = ", 
  round(sqrt(sum((sim_data$dist - dist(sim_res$proposed$lambda0.7$z))^2) /
               sum((dist(sim_res$proposed$lambda0.7$z))^2)), 4)), side=3)
dev.off()


## Shepard
pdf("data/result/shepard_sim.pdf",
    width = 9, height = 9)
par(mfrow = c(2, 2))
plot(sim_data$dist, dist(sim_res$mds),
     xlab = "original distance", ylab = "configuration distance",
     main = "Pure MDS")
mtext(paste(
  "Stress1 = ", 
  round(sqrt(sum((sim_data$dist - dist(sim_res$mds))^2) /
               sum((dist(sim_res$mds))^2)), 4)), side=3)

plot(sim_data$dist, dist(sim_res$proposed$lambda0.3$z),
     xlab = "original distance", ylab = "configuration distance",
     main = expression(paste(lambda, "= 0.3")))
mtext(paste(
  "Stress1 = ", 
  round(sqrt(sum((sim_data$dist - dist(sim_res$proposed$lambda0.3$z))^2) /
               sum((dist(sim_res$proposed$lambda0.3$z))^2)), 4)), side=3)

plot(sim_data$dist, dist(sim_res$proposed$lambda0.5$z),
     xlab = "original distance", ylab = "configuration distance",
     main = expression(paste(lambda, "= 0.5")))
mtext(paste(
  "Stress1 = ", 
  round(sqrt(sum((sim_data$dist - dist(sim_res$proposed$lambda0.5$z))^2) /
               sum((dist(sim_res$proposed$lambda0.5$z))^2)), 4)), side=3)

plot(sim_data$dist, dist(sim_res$proposed$lambda0.7$z),
     xlab = "original distance", ylab = "configuration distance",
     main = expression(paste(lambda, "= 0.7")))
mtext(paste(
  "Stress1 = ", 
  round(sqrt(sum((sim_data$dist - dist(sim_res$proposed$lambda0.7$z))^2) /
               sum((dist(sim_res$proposed$lambda0.7$z))^2)), 4)), side=3)
dev.off()

