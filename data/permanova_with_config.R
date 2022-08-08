source("get_dist.R")

ordu1 = ordinate(site1, "PCoA", distance = "unifrac", weighted=TRUE)
ordu2 = ordinate(site2, "PCoA", distance = "unifrac", weighted=TRUE)

pseudo_F <- function(mat=NULL, trt, d = NULL){
  if(is.null(d)){
    d <- as.matrix(dist(mat))
  } else {
    d <- as.matrix(d)
  }
  N <- nrow(d)
  a <- length(unique(trt))
  n <- N/a
  
  SST <- SSW <- 0
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      SST = SST + d[i,j]^2
      SSW <- SSW + (trt[i] == trt[j]) * d[i,j]^2
    }
  }
  SST <- SST / N
  SSW <- SSW / n
  SSA <- SST - SSW
  return(list(SST = SST, SSW = SSW, SSA = SSA,
              pseudoF = (SSA * (N-a))/(SSW * (a-1))))
}

pseudo_F(ordu1$vectors[,1:2], site1@sam_data@.Data[[1]])
pseudo_F(ordu2$vectors[,1:2], site2@sam_data@.Data[[1]])