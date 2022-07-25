library("phyloseq")
library("ggplot2")
library("ape")
library("vegan")
library(dplyr)

#### Load Data
ps = readRDS("community_phyloseq.Rds")


#### Distance matrix
site1 = subset_samples(ps, Site=="1") # subset
site2 = subset_samples(ps, Site=="2") # subset
dist1 = phyloseq::distance(site1, method = "unifrac", weighted = T)
dist2 = phyloseq::distance(site2, method = "unifrac", weighted = T)

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


#### Plot & Save Distance matrix
pdf("dist_heatmap.pdf")
data.frame(x = rep(1:36, each = 36),
           y = rep(1:36, 36),
           val = c(distmat1)) %>%
  ggplot(., aes(x=x, y=y, fill = val)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="blue") +
  geom_vline(aes(xintercept = 18.5), colour = "red", size = 1) +
  geom_hline(aes(yintercept = 18.5), color = "red", size = 1) +
  annotate("text", label = c("Pt+", "Pt-", "Pt+", "Pt-"),
           x = c(10, 28, 0, 0), y = c(0,0, 10, 28), color = "red", size = 6) +
  ggtitle("Site 1")

data.frame(x = rep(1:36, each = 36),
           y = rep(1:36, 36),
           val = c(distmat2)) %>%
  ggplot(., aes(x=x, y=y, fill = val)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="blue") +
  geom_vline(aes(xintercept = 18.5), colour = "red", size = 1) +
  geom_hline(aes(yintercept = 18.5), color = "red", size = 1) +
  annotate("text", label = c("Pt+", "Pt-", "Pt+", "Pt-"),
           x = c(10, 28, 0, 0), y = c(0,0, 10, 28), color = "red", size = 6) +
  ggtitle("Site 2")
dev.off()
