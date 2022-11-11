## Created and last run 2/13/22 by Hyu Kim (hskimm@mit.edu)
library("phyloseq")
library("ggplot2")
library("ape")
library("vegan")

# path = "/Users/hk/desktop/mit/doctoral/idps/2021-08" # write your current directory here
ps = readRDS("community_phyloseq.Rds")

table(ps@sam_data$Site)
# ps@tax_table
table(ps@sam_data$Site, ps@sam_data$Treatment)


###### plot PCoA ######
# "Unifrac" distance metric, see: https://doi.org/10.1128/AEM.71.12.8228-8235.2005
ps_filter = subset_samples(ps, Site=="1") # subset each site ("1" or "2")
ordu = ordinate(ps_filter, "PCoA", distance = "unifrac", weighted=TRUE);

sum(ordu$values$Eigenvalues > 0)

#in ordu, only e-vectors of e-val>0 are saved.
#to get original data, need to see unifrac(metric matrix: not Sigma as usual PCA.)


# zz <- plot_ordination(ps_filter, ordu, shape="Treatment", axes = 1:10)
# ord_plot$data$Axis.1
# zz$data$Axis.1
# ordu$values$Eigenvalues
# dim(ordu$vectors)
# 
# ordu$values$Cum_corr_eig
# 
# library("plot3D")
# 
# scatter3D(zz$data$Axis.1,
#           zz$data$Axis.2,
#           zz$data$Axis.3,
#           pch = 19, col = ifelse(zz$data$Treatment == "Pt +", "red", "black"))
# 
# zz$data$Treatment
# 
# # install.packages("devtools")  # so we can install from github
# library("devtools")
# # install_github("ropensci/plotly")  # plotly is part of ropensci
# library(plotly)
# fig <- plot_ly(zz$data,
#                x = ~Axis.1, y = ~Axis.2, z = ~Axis.3,
#                color = ~Treatment, colors = c("red", "black"))
# fig
# 
# ps_filter2 = subset_samples(ps, Site=="2") # subset each site ("1" or "2")
# ordu2 = ordinate(ps_filter2, "PCoA", distance = "unifrac", weighted=TRUE);
# zz2 <- plot_ordination(ps_filter2, ordu2, shape="Treatment", axes = 1:10)
# fig2 <- plot_ly(zz2$data,
#                x = ~Axis.1, y = ~Axis.2, z = ~Axis.3,
#                color = ~Treatment, colors = c("red", "black"))
# fig2


ord_plot <- plot_ordination(ps_filter, ordu, shape="Treatment", axes = 1:10) +
  geom_point(aes(shape=Treatment), size=1.35, alpha=2) + 
  scale_shape_manual(values = c(1,16)) + # label condition color
  # scale_y_continuous(breaks=c(-0.05, 0, 0.05)) + # for Site 1 only
  theme(text = element_text(size=8), 
        legend.position = "none",
        # axis.text.y   = element_blank(),
        # axis.text.x   = element_blank(),
        # axis.title.y  = element_blank(),
        # axis.title.x  = element_blank(),
        axis.text.y   = element_text(size=8, colour="black"),
        axis.text.x   = element_text(size=8, colour="black"),
        axis.title.y  = element_text(size=7, colour="black"),
        axis.title.x  = element_text(size=7, colour="black"),
        strip.background = element_rect(fill=NA),
        strip.text.x = element_blank(),
        axis.ticks = element_line(colour = "black", size = 0.25),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        # axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.3))

plot(ord_plot)

# setwd("/Users/hk/Desktop/MIT/Doctoral/IDPS")
# ggsave("plot3.eps", width = 1.8, height = 1.7, units = "in")


###### Statistical test "PERMANOVA" ######
# Perform test between two groups, with / without host, on each site ("Layer")
# Based on Bray–Curtis dissimilarity or Unifrac distance metric
for (j in c("1","2")){
  physeq_temp = subset_samples(ps, Site==j) # subset
  bray = phyloseq::distance(physeq_temp, method = "unifrac")
  sampledf = data.frame(sample_data(physeq_temp))
  
  ## 1) Adonis test (compare centroids between each group)
  ado = adonis(bray ~ Treatment, data = sampledf, strata=sampledf$Replicate)
  cat("\n$$$$$ Site ", j, " $$$$$")
  print(ado)
  
  ## 2) Homogeneity of dispersion test (compare dispersion of the groups)
  beta <- betadisper(bray, sampledf$Treatment)
  print(permutest(beta))
}

site1 = subset_samples(ps, Site=="1") # subset
site2 = subset_samples(ps, Site=="2") # subset

site1@sam_data$Treatment
site2@sam_data$Treatment


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



library(dplyr)

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
