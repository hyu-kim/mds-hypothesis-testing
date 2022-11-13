# Objects dependency on 'rigorousness_check.R'
# Within a set of randomly labeled ps dataset, finds cases
# where P value is affected by MDS reduction

source("permanova_with_config.R")
source('rigorousness_check.R')


# function
pick_extreme <- function(random_f, labels){
  ind_significant <- (random_f[,3]>2)|(random_f[,4]>2)
  random_f_filtered <- random_f[ind_significant,]
  p_ratio <- random_f_filtered[,3]/random_f_filtered[,4]
  index_extremes <- which(p_ratio==max(p_ratio) | p_ratio==min(p_ratio))
  labels_extremes <- labels[ind_significant,][index_extremes,]
  
  return(list(f=random_f_filtered, index=index_extremes, label=labels_extremes))
}


# call function
extremes1_list <- pick_extreme(random_f1, labels1)
extremes2_list <- pick_extreme(random_f2, labels2)
rand_f1_sub <- extremes1_list$f
rand_f2_sub <- extremes2_list$f
index_extreme1 <- extremes1_list$index
index_extreme2 <- extremes2_list$index
label_extreme1 <- extremes1_list$label
label_extreme2 <- extremes2_list$label

zmds1 <- plot_ordination(ps1, ordu1, axes = 1:10)$data[, 1:2]
zmds2 <- plot_ordination(ps2, ordu2, axes = 1:10)$data[, 1:2]

# plot MDS
par(mfrow = c(1,2))
plot(zmds1, col = label_extreme1[1,], 
     main = paste("P0=", 10^-rand_f1_sub[index_extreme1[1],3],
                  "P_mds=", 10^-rand_f1_sub[index_extreme1[1],4]))
plot(zmds1, col = label_extreme1[2,], 
     main = paste("P0=", 10^-rand_f1_sub[index_extreme1[2],3],
                  "P_mds=", 10^-rand_f1_sub[index_extreme1[2],4]))

par(mfrow = c(1,2))
plot(zmds2, col = label_extreme2[1,], 
     main = paste("P0=", 10^-rand_f2_sub[index_extreme2[1],3],
                  "P_mds=", 10^-rand_f2_sub[index_extreme2[1],4]))
plot(zmds2, col = label_extreme2[2,], 
     main = paste("P0=", 10^-rand_f2_sub[index_extreme2[2],3],
                  "P_mds=", 10^-rand_f2_sub[index_extreme2[2],4]))