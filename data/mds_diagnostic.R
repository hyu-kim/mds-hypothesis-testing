# variables dependency on get_dist.R
source("mm.R")

#### Shepard Plot
# Some definitions switches x and y axis
par(mfrow = c(1,2))
plot(dist2, dist(zmds2),
     xlab = "original distance", ylab = "configuration distance",
     main = "Shepard Diagram Pure MDS")
mtext(paste(
  "Stress1 = ", 
  round(sqrt(sum((dist2 - dist(obmm_x$z))^2) /
                   sum((dist(obmm_x$z))^2)), 4)), side=3)
plot(dist2, dist(obmm_x$z),
     xlab = "original distance", ylab = "configuration distance",
     main = expression(paste("Shepard Diagram ", lambda, "= 0.3")))
mtext(paste(
  "Stress1 = ", 
  round(sqrt(sum((dist2 - dist(zmds2))^2) /
                   sum((dist(zmds2))^2)), 4)), side=3)
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

ggplot(data = data.frame(obmm_x$z, y = as.factor(y2s[,1]), spp = spp()),
       aes(x = Axis.1, y = Axis.2, size = spp,
           color = y)) +
  geom_point() +
  labs(title="Stress-Per-Point Bubble Plot")
