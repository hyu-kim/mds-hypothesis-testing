library("phyloseq")

ps = readRDS("community_phyloseq.Rds")

df <- psmelt(ps)

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

df_summ <- add_taxa(df)

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

df_out <- get_matrix(df_summ)
write.table(df_out, file = "result/Alga/abundance.tsv", row.names=TRUE, col.names=FALSE, sep="\t")

lab <- as.numeric(grepl("Con1", colnames(df_out)))
write.table(lab, file = "result/Alga/labels.txt", sep='\n', col.names=FALSE, row.names=FALSE)