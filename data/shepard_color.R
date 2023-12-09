# source("mm.R")
library(superMDS)
library(scales)
library(grid)
library(gridExtra)
library(parallel)
library("ggplot2")

load('result/sim_etc.RData')

# blue 0944BB
# oran DB4915
# lightblue BCDDF5
# ligntoran F8CEBF

## Simulated superMDS
print(Sys.time())
sim_res_smds <- mclapply(1:4, function(i){
  x <- c(0, 0.3, 0.5, 0.7)[i]
  return(TrainSuperMDS(d = sim_data$distmat, y = sim_data$data$Y, alpha = x/(1+x)))
  print(Sys.time())
}, mc.cores = 4)
names(sim_res_smds) <- c("lambda0", "lambda0.3", "lambda0.5", "lambda0.7")


## Shepard with color
pdf("result/shepard_sim_colored.pdf", width = 7, height = 8)
par(mfrow = c(2, 2))

plot(sim_data$dist, dist(sim_res$mds),
     xlab = "original", ylab = "configuration",
     xaxt="n", yaxt="n",
     main = "proposed, lambda = 0",
     col = alpha("#0944BB", 0.2), 
     pch=16
     )
axis(side=1, at=seq(0,10,5), labels = FALSE)
axis(side=2, at=seq(0,10,5), labels = FALSE)

plot(sim_data$dist, dist(sim_res$proposed$lambda0.7$z),
     xlab = "original", ylab = "configuration",
     xaxt="n", yaxt="n",
     main = "proposed, lambda = 0.7", 
     col = alpha("#0944BB", 0.2), 
     pch=16
     )
axis(side=1, at=seq(0,10,5), labels = FALSE)
axis(side=2, at=seq(0,10,5), labels = FALSE)

plot(sim_data$dist, dist(sim_res_smds$lambda0$z),
     xlab = "original", ylab = "configuration",
     xaxt="n", yaxt="n",
     main = "SMDS, lambda = 0", 
     col = alpha("#DB4915", 0.2), 
     pch=16
     )
axis(side=1, at=seq(0,10,5), labels = FALSE)
axis(side=2, at=seq(0,10,5), labels = FALSE)

plot(sim_data$dist, dist(sim_res_smds$lambda0.7$z),
     xlab = "original", ylab = "configuration",
     xaxt="n", yaxt="n",
     main = "SMDS, lambda = 0.7", 
     col = alpha("#DB4915", 0.15), 
     pch=16
     )
axis(side=1, at=seq(0,10,5), labels = FALSE)
axis(side=2, at=seq(0,10,5), labels = FALSE)

dev.off()


## Configuration - brought from pcoa.R
ord_plot <- ggplot(res1$lambda0.5$z, aes(Axis.1, Axis.2)) + # choose one either below
# ord_plot <- ggplot(zmds2, aes(Axis.1, Axis.2)) +
  geom_point(aes(shape=as.factor(y1)), size=2.5, stroke=0.5, colour="#DB4915", fill="#F8CEBF") +
  scale_shape_manual(values = c(16,21)) +
  # scale_y_continuous(breaks=c(-0.05, 0, 0.05)) + # Layer 1 only
  theme(text = element_text(size=7), 
        legend.position = "none",
        axis.title.y  = element_blank(),
        axis.title.x  = element_blank(),
        axis.text.y   = element_blank(),
        axis.text.x   = element_blank(),
        strip.background = element_rect(fill=NA),
        strip.text.x = element_blank(),
        axis.ticks = element_line(colour = "black", size = 0.25),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        # axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))

setwd("~/Desktop/MIT/Doctoral/defense/materials/codes")
ggsave("pcoa_site1_prop.eps", width = 2, height = 2, units = "in")