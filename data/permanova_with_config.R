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
              ratio = (SSA * (N-a))/(SSW * (a-1))))
}


get_p <- function(mat=NULL, d=NULL, trt, n_iter=999, fun=pseudo_F){
  # initialize
  f_permuted = matrix(0, nrow=n_iter, ncol=1)  # pseudo-F only
  # iterate to get pseudo F
  N <- length(trt)
  tbl <- table(trt)
  for (iter in 1:n_iter){
    ind_rand <- sample(1:N, tbl[2], replace=F)
    y_rand <- rep(1,N)
    y_rand[ind_rand] = 2
    f_permuted[iter,1] = fun(mat=mat, d = d, trt = y_rand)$ratio
  }
  # compute p value
  f_sorted = sort(f_permuted[,1], decreasing = TRUE)
  f_val = fun(mat=mat, d = d, trt = trt)$ratio
  p_val = which(f_val > f_sorted)[1]
  p_val <- (p_val-1)/(n_iter+1)
  
  return(list(ratio_all = f_sorted, ratio = f_val, p = p_val))
}

pseudo_F(ordu1$vectors[,1:2], site1@sam_data@.Data[[1]])
pseudo_F(ordu2$vectors[,1:2], site2@sam_data@.Data[[1]])