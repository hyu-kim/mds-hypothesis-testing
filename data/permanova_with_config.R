library(sys)
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


pseudo_F_iter <- function(mat=NULL, d=NULL, y_iter){
  # mat   : Configuration. Typically N by 2 matrix
  # d     : distance matrix, square, N by N
  # y_iter: iterated labels, N by n_iter
  N <- nrow(y_iter)
  n_iter <- ncol(y_iter)
  
  if(is.null(d)){
    d <- as.matrix(dist(mat))
  } else {
    d <- as.matrix(d)
  }
  
  # create an iterated indicator matrix, sized N by N by n_iter
  # ind_iter <- array(0, dim=c(N, N, n_iter))
  # for (it in 1:n_iter){
  y_iter_3d <- array(rep(y_iter, N), dim=c(N, n_iter, N))
  y_iter_3d_t <- aperm(y_iter_3d, c(3,2,1))
  ind_iter <- (y_iter_3d == y_iter_3d_t) * 1
  ind_iter <- aperm(ind_iter, c(1,3,2))  # now (N, N, n_iter)
  # }
  d_3d <- array(rep(d, n_iter), dim=c(N, N, n_iter))
  nomi <- 2 * rowSums(aperm(ind_iter * d_3d * d_3d, c(3,2,1)), dims=1)  # merge over N by N, leaving n_iter
  denomi <- d_sq_sum <- sum(d^2)  # this is denominator
  return((denomi/nomi - 1)*(N-2))
}


get_p <- function(mat=NULL, d=dist1, trt=y1, n_iter=999, fun=pseudo_F){
  time1 = Sys.time()
  
  # initialize then iterate to get pseudo F
  # f_permuted = matrix(0, nrow=n_iter, ncol=1)  # old version
  N <- length(trt)
  tbl <- table(trt)
  ind_iter <- replicate(n=n_iter, sample(1:N, tbl[2], replace=F))  # new version
  y_iter <- matrix(1, nrow=N, ncol=n_iter)  # new version
  for (iter in 1:n_iter){
    y_iter[ind_iter[,iter],iter] <- 2  # new version
    # ind_rand <- sample(1:N, tbl[2], replace=F)  # old version
    # y_rand <- rep(1,N)  # old version
    # y_rand[ind_rand] = 2  # old version
    # f_permuted[iter,1] = fun(mat=mat, d = d, trt = y_rand)$ratio  # old version
  }
  f_permuted <- pseudo_F_iter(mat=mat, d=d, y_iter=y_iter)  # new version
  # compute p value
  f_sorted = sort(f_permuted, decreasing = TRUE)
  f_val = fun(mat=mat, d = d, trt = trt)$ratio
  p_val = which(f_val > f_sorted)[1]
  p_val <- (p_val-1)/(n_iter+1)
  
  time2 = Sys.time()
  # print(paste('elapsed:', time2-time1))
  return(list(ratio_all = f_sorted, ratio = f_val, p = p_val))
}

pseudo_F(ordu1$vectors[,1:2], site1@sam_data@.Data[[1]])
pseudo_F(ordu2$vectors[,1:2], site2@sam_data@.Data[[1]])