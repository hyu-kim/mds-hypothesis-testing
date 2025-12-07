# loads sim_1_rev data
# permute labels, compute F ratios in Z and X
# compare F ratios ordered vs. as-is, two panels
# compare F ratios permuted once for all vs. twice for each
# highlight point where un-permuted

library(vegan)
library(scales)
source("permanova_with_config.R")

get_y_rand <- function(y){
  N <- 200
  a <- 2
  tbl <- table(y)
  y_rand <- rep(1,N)
  pool_v <- 1:N
  for(cl in 2:a){
    ind_rand <- sample(pool_v, tbl[cl], replace=F)
    y_rand[ind_rand] = cl
    pool_v <- setdiff(pool_v, ind_rand)
  }
  
  return(y_rand)
}

get_f_pair <- function(x_mat, y){
  d_x <- vegdist(t(x_mat), method="bray")
  pcoa <- cmdscale(d_x, eig = TRUE)
  pcoa_df <- as.data.frame(pcoa$points)
  d_z <- dist(pcoa_df, method='euclidean')
  f_x <- pseudo_F(d = d_x, trt = y)
  f_z <- pseudo_F(d = d_z, trt = y)
  
  return(list(f_x = f_x, f_z = f_z))
}

permute_pair_f <- function(n_iter, x_mat, y){
  f_paired_mat = matrix(0, nrow=n_iter+1, ncol=2) # 1st, x; 2nd, z
  f_paired <- get_f_pair(x_mat, y)
  f_paired_mat[1,1]<- f_paired$f_x
  f_paired_mat[2,2]<- f_paired$f_z
  
  for(i in 2:(n_iter+1)){
    y_rand <- get_y_rand(y)
    f_paired <- get_f_pair(x_mat, y_rand)
    f_paired_mat[i,1]<- f_paired$f_x
    f_paired_mat[i,2]<- f_paired$f_z
  }
  
  return(f_paired_mat)
}

permute_once_sort_f <- function(n_iter, x_mat, y){ # once for both x and z
  f_paired_mat <- permute_pair_f(n_iter, x_mat, y)  
  f_sorted_x = sort(f_paired_mat[2:(n_iter+1),1], decreasing = TRUE)
  f_sorted_z = sort(f_paired_mat[2:(n_iter+1),2], decreasing = TRUE)
  f_sorted_mat <- cbind(f_sorted_x, f_sorted_z)
  
  return(f_sorted_mat)
}

permute_each_sort_f <- function(n_iter, x_mat, y){ # independently for x and z
  f_sorted_mat_each = matrix(0, nrow=n_iter, ncol=2)
  set.seed(100)
  f_sorted_mat_each[,1] <- permute_once_sort_f(n_iter, x_mat, y)[,1]
  set.seed(101)
  f_sorted_mat_each[,2] <- permute_once_sort_f(n_iter, x_mat, y)[,2]
  
  return(f_sorted_mat_each)
}

x_mat <- readRDS('result/ScalingStudy/sim_rev_1/sim_rev_1-N200-data.Rds')
y <- read.csv('result/ScalingStudy/sim_rev-N200-Y.csv')$x

n_iter <- 999

par(mfrow = c(1, 3))

f_paired_mat <- permute_pair_f(n_iter, x_mat, y)
plot(f_paired_mat[,1], f_paired_mat[,2], xlab="F_x", ylab="F_z", pch=16, col = alpha("black", 0.5))
points(f_paired_mat[1,1], f_paired_mat[1,2], col = "red", pch = 4, cex = 2)

f_sorted_mat <- permute_once_sort_f(n_iter, x_mat, y)
plot(f_sorted_mat[,1], f_sorted_mat[,2])
points(f_paired_mat[1,1], f_paired_mat[1,2], col = "red", pch = 4)

f_each_sorted_mat <- permute_each_sort_f(n_iter, x_mat, y)
plot(f_each_sorted_mat[,1], f_each_sorted_mat[,2])
points(f_paired_mat[1,1], f_paired_mat[1,2], col = "red", pch = 4)