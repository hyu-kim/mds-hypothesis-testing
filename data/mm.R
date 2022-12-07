# variables dependency on get_dist.R, data.R
source("permanova_with_config.R")


# Distance between vector
get_dist_mat <- function(z){
  N = dim(z)[1]
  z_dist = matrix(0, nrow = N, ncol = N)
  for(i in 1:N){
    for(j in 1:N){
      z_dist[i,j] <- sqrt(sum((z[i,] - z[j,])^2))
    }
  }
  return(z_dist)
}


# Indicator matrix of labels
get_ind_mat <- function(y){
  N = length(y)
  mat = matrix(0, nrow = N, ncol = N)
  for(i in 1:N){
    for(j in 1:N){
      mat[i,j] <- as.numeric(y[i]==y[j])
    }
  }
  return(mat)
}


# MDS objective
mds_obj <- function(D, z){
  z_distmat <- get_dist_mat(z)
  return(sum((D - z_distmat)^2)/2)
}


# confirmatory objective term with labels
conf_obj <- function(y, z, D){
  z_distmat <- get_dist_mat(z)
  y_indmat <- get_ind_mat(y)
  phi <- sum((1-y_indmat) * D*D) / sum(y_indmat * D*D)
  val <- (1-(1+phi)*y_indmat) * z_distmat^2
  res <- 0.5 * abs(sum(val))
  return(list(val = res, sign = sign(rowSums(val))))
}


mm_cmds <- function(nit = 100, conv_crit = 5e-03, lambda = 0.2,
                    z0, D, y){
  N <- dim(z0)[1]
  S <- dim(z0)[2]
  y_indmat <- get_ind_mat(y)
  phi <- sum((1-y_indmat) * D*D) / sum(y_indmat * D*D)
  delta <- conf_obj(y, z0, D)$sign
  z_temp <- z_up <- z0
  for(t in 1:nit){
    obj_conf_up <- conf_obj(y, z_up, D)
    obj_up <- obj_conf_up$val + lambda*mds_obj(D, z_up)
    print(paste('iteration', t, 
                '  total', sprintf(obj_up, fmt = '%#.3f'), 
                '  conf_term', sprintf(obj_conf_up$val, fmt = '%#.3f')
    ))
    
    delta <- obj_conf_up$sign  # N
    z_distmat <- get_dist_mat(z_up)  # (N,N)
    coeff <- D/z_distmat  # final term in the update
    coeff[is.nan(coeff)] <- 0
    for(i in 1:N){
      z_diff <- -sweep(x=z_up, MARGIN=2, STATS=as.matrix(z_up[i,]), FUN="-")
      z_temp[i,] <- (1+lambda*delta[i])*apply(z_up[y!=y[i],], 2, sum) +
        (1-lambda*phi*delta[i])*apply(z_up[y==y[i],], 2, sum) +
        sweep(x=z_diff, MARGIN=1, STATS=coeff[,i], FUN="*")
      z_temp[i,] <- z_temp[i,] / (2 + lambda * (1-phi)*delta[i])
    }
    z_up <- z_temp
  }
  return(list(z = z_up, obj0 = obj0, obj_up = obj_up,
              F0 = F0, Finit = Finit, Fup = Fz_cur))
}