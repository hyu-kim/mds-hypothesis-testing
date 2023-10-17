## creates and exports dataframe for neural network training
library("phyloseq")

ps = readRDS("community_phyloseq.Rds")
site1 = subset_samples(ps, Site=="1") # subset
site2 = subset_samples(ps, Site=="2")

site1_df <- psmelt(site1)
site2_df <- psmelt(site2)

add_taxa <- function(df, tax_abb=c('k','p','c','o','f','g','s')){
  df_summ <- df[,c("Sample",'Abundance')]
  df_summ$Taxon <- NA
  for (i in seq(nrow(df))){
    taxon <- ''
    for (j in 1:6){
      taxon <- paste(taxon, '|', tax_abb[j], '__', df[i,j+7], sep='')
    }
    taxon <- paste(taxon, '|', tax_abb[7], '__', df$OTU[i], sep='')
    taxon <- substr(taxon, 2, nchar(taxon))
    df_summ[i,'Taxon'] <- taxon
  }
  return(df_summ)
}

site1_df_summ <- add_taxa(site1_df)
site2_df_summ <- add_taxa(site2_df)

get_matrix <- function(df_summ){
  df_out <- data.frame(matrix(0, nrow=length(unique(df_summ$Taxon)), ncol=length(unique(df_summ$Sample))))
  rownames(df_out) <- unique(df_summ$Taxon)
  colnames(df_out) <- unique(df_summ$Sample)
  for (i in seq(nrow(df_out))){  # per taxon
    ta <- rownames(df_out)[i]
    df_summ_f <- df_summ[df_summ$Taxon==ta,]
    for (j in seq(nrow(df_summ_f))){  # per sample
      sa <- unique(df_summ_f$Sample)[j]
      ind <- which(colnames(df_out)==sa)
      df_out[i,ind] <- df_summ_f[j, 'Abundance']
    }
  }
  return(df_out)
}

site1_df_out <- get_matrix(site1_df_summ)
site2_df_out <- get_matrix(site2_df_summ)

setwd('..')

write.table(site1_df_out, file = "PopPhy/data/Alga/abundance_site1.tsv", 
            row.names=TRUE, col.names=FALSE, sep="\t")
write.table(site2_df_out, file = "PopPhy/data/Alga/abundance_site2.tsv", 
            row.names=TRUE, col.names=FALSE, sep="\t")

lab1 <- as.numeric(grepl("Con1", colnames(site1_df_out)))
lab2 <- as.numeric(grepl("Con1", colnames(site2_df_out)))


write.table(lab1, file = "PopPhy/data/Alga/labels_site1.txt", sep='\n', 
            col.names=FALSE, row.names=FALSE)
write.table(lab2, file = "PopPhy/data/Alga/labels_site2.txt", sep='\n', 
            col.names=FALSE, row.names=FALSE)