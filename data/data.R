library("phyloseq")
library("ggplot2")
library("ape")
library("vegan")

# path = "/Users/hk/desktop/mit/doctoral/idps/2021-08" # write your current directory here
ps = readRDS("community_phyloseq.Rds")

ps_filter1 = subset_samples(ps, Site=="1") # subset each site ("1" or "2")
ordu1 = ordinate(ps_filter1, "PCoA", distance = "unifrac", weighted=TRUE);
zz1 <- plot_ordination(ps_filter1, ordu1, shape="Treatment", axes = 1:10)
zmds1 <- zz1$data[, 1:2]

ps_filter2 = subset_samples(ps, Site=="2") # subset each site ("1" or "2")
ordu2 = ordinate(ps_filter2, "PCoA", distance = "unifrac", weighted=TRUE);
zz2 <- plot_ordination(ps_filter2, ordu2, shape="Treatment", axes = 1:10)
zmds2 <- zz2$data[, 1:2]