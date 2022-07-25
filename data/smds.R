source("get_dist.R")
library(superMDS)

y1 <- ifelse(site1@sam_data$Treatment == "Pt +", 1, 2)
y2 <- ifelse(site2@sam_data$Treatment == "Pt +", 1, 2)

#### SMDS to our data
smds_site1 <- list(
  a.0 =  TrainSuperMDS(d = distmat1, y = y1, alpha = 0.0),
  a.1 =  TrainSuperMDS(d = distmat1, y = y1, alpha = 0.1),
  a.2 =  TrainSuperMDS(d = distmat1, y = y1, alpha = 0.2),
  a.3 =  TrainSuperMDS(d = distmat1, y = y1, alpha = 0.3),
  a.4 =  TrainSuperMDS(d = distmat1, y = y1, alpha = 0.4),
  a.5 =  TrainSuperMDS(d = distmat1, y = y1, alpha = 0.5),
  a.6 =  TrainSuperMDS(d = distmat1, y = y1, alpha = 0.6),
  a.7 =  TrainSuperMDS(d = distmat1, y = y1, alpha = 0.7),
  a.8 =  TrainSuperMDS(d = distmat1, y = y1, alpha = 0.8),
  a.9 =  TrainSuperMDS(d = distmat1, y = y1, alpha = 0.9),
  a.10 = TrainSuperMDS(d = distmat1, y = y1, alpha = 1)
)
smds_site2 <- list(
  a.0 =  TrainSuperMDS(d = distmat2, y = y2, alpha = 0.0),
  a.1 =  TrainSuperMDS(d = distmat2, y = y2, alpha = 0.1),
  a.2 =  TrainSuperMDS(d = distmat2, y = y2, alpha = 0.2),
  a.3 =  TrainSuperMDS(d = distmat2, y = y2, alpha = 0.3),
  a.4 =  TrainSuperMDS(d = distmat2, y = y2, alpha = 0.4),
  a.5 =  TrainSuperMDS(d = distmat2, y = y2, alpha = 0.5),
  a.6 =  TrainSuperMDS(d = distmat2, y = y2, alpha = 0.6),
  a.7 =  TrainSuperMDS(d = distmat2, y = y2, alpha = 0.7),
  a.8 =  TrainSuperMDS(d = distmat2, y = y2, alpha = 0.8),
  a.9 =  TrainSuperMDS(d = distmat2, y = y2, alpha = 0.9),
  a.10 = TrainSuperMDS(d = distmat2, y = y2, alpha = 1)
)

sapply(1:10, function(x){min(smds_site1[[x]]$crits)})
sapply(1:10, function(x){min(smds_site2[[x]]$crits)})


#### Save plots
pdf("SMDS_res.pdf", width = 7.5, height = 5)
par(mfrow = c(2, 3), mar = c(2,2,2,1))
plot(smds_site1$a.0$z, pch = ifelse(y1 == 1, 16, 1), cex = 1.5,
     xlab = "", ylab = "", main = "Site 1, alpha = 0")
plot(smds_site1$a.5$z, pch = ifelse(y1 == 1, 16, 1), cex = 1.5,
     xlab = "", ylab = "", main = "Site 1, alpha = 0.5")
plot(smds_site1$a.10$z, pch = ifelse(y1 == 1, 16, 1), cex = 1.5,
     xlab = "", ylab = "", main = "Site 1, alpha = 1")
legend("bottomright", c("Pt +", "Pt -"), pch = c(16, 1))
plot(smds_site2$a.0$z, pch = ifelse(y1 == 1, 16, 1), cex = 1.5,
     xlab = "", ylab = "", main = "Site 2, alpha = 0")
plot(smds_site2$a.5$z, pch = ifelse(y1 == 1, 16, 1), cex = 1.5,
     xlab = "", ylab = "", main = "Site 2, alpha = 0.5")
plot(smds_site2$a.10$z, pch = ifelse(y1 == 1, 16, 1), cex = 1.5,
     xlab = "", ylab = "", main = "Site 2, alpha = 1")
legend("bottomright", c("Pt +", "Pt -"), pch = c(16, 1))
dev.off()


pdf("SMDS_crit.pdf", width = 7.5, height = 2.5)
par(mfrow = c(1,3), mar = c(3,3,2,1), mgp=c(2,1,0))
plot((0:10)/10, 
     sapply(1:11, function(x){min(smds_site1[[x]]$crits)}),
     ylim = c(0.05, 0.6),
     type = "o", col = "blue",
     xlab = "alpha", ylab = "criterion",
     main = "Minimized Criterion")
points((0:10)/10, 
     sapply(1:11, function(x){min(smds_site2[[x]]$crits)}),
     type = "o", col = "red", lty = 2)

plot((0:10)/10, 
     sapply(1:11, function(x){smds_site1[[x]]$stress}),
     type = "o", col = "blue",
     ylim = c(0.2, 1.5),
     xlab = "alpha", ylab = "stress",
     main = "Stress function")
points((0:10)/10, 
       sapply(1:11, function(x){smds_site2[[x]]$stress}),
       type = "o", col = "red", lty = 2)

plot((0:10)/10, 
     sapply(1:11, function(x){smds_site1[[x]]$super}),
     ylim = c(0, 5.5),
     type = "o", col = "blue",
     xlab = "alpha", ylab = "supervised",
     main = "Supervised term")
points((0:10)/10, 
       sapply(1:11, function(x){smds_site2[[x]]$super}),
       type = "o", col = "red", lty = 2)
legend("topright", c("Site 1", "Site 2"), 
       col = c("blue", "red"),
       lty = c(1,2))
dev.off()
