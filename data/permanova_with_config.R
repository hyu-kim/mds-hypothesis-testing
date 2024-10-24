# source("get_dist.R")

# ordu1 = ordinate(site1, "PCoA", distance = "unifrac", weighted=TRUE)
# ordu2 = ordinate(site2, "PCoA", distance = "unifrac", weighted=TRUE)
get_ind_mat <- function(y){
  N = length(y)
  mat = matrix(0, nrow = N, ncol = N)
  for(i in 1:N){
    for(j in 1:(i-1)){
      mat[i,j] <- as.numeric(y[i]==y[j])
    }
  }
  return(mat + t(mat))
}


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
  SST <- sum(d * d) / 2
  SSW <- sum(d * d * get_ind_mat(trt)) / 2
  SST <- SST / N
  SSW <- SSW / n
  SSA <- SST - SSW
  return((SSA * (N-a))/(SSW * (a-1)))
}


get_p <- function(mat=NULL, d=NULL, trt, n_iter=999, fun=pseudo_F){
  # initialize
  f_permuted = matrix(0, nrow=n_iter, ncol=1)  # pseudo-F only
  # iterate to get pseudo F
  N <- length(trt)
  a <- length(unique(trt))
  tbl <- table(trt)
  for (iter in 1:n_iter){
    y_rand <- rep(1,N)
    pool_v <- 1:N
    for(cl in 2:a){
      ind_rand <- sample(pool_v, tbl[cl], replace=F)
      y_rand[ind_rand] = cl
      pool_v <- setdiff(pool_v, ind_rand)
    }
    f_permuted[iter,1] = fun(mat=mat, d = d, trt = y_rand)
  }
  # compute p value
  f_sorted = sort(f_permuted[,1], decreasing = TRUE)
  f_val = fun(mat=mat, d = d, trt = trt)
  p_val = which(f_val > f_sorted)[1]
  p_val <- (p_val-1)/(n_iter+1)
  
  return(list(ratio_all = f_sorted, ratio = f_val, p = p_val))
}

# pseudo_F(ordu1$vectors[,1:2], site1@sam_data@.Data[[1]])
# pseudo_F(ordu2$vectors[,1:2], site2@sam_data@.Data[[1]])