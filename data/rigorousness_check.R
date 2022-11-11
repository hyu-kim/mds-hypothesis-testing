# Randomly assigns a binary label to a given phyloseq dataset theb calculates
# pseudo-F under original and MDS ordination respectively, which is repeated
# for a number of times
 
library("phyloseq")
source("permanova_with_config.R")

ps = readRDS("community_phyloseq.Rds")

