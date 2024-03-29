# variables dependency on get_dist.R
source("permanova_with_config.R")

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


# compute phi (ratio)
get_phi <- function(mat=NULL, d=NULL, trt){
  if(is.null(d)){
    d <- as.matrix(dist(mat))
  } else {
    d <- as.matrix(d)
  }
  y_indmat <- get_ind_mat(trt)
  phi <- sum((1-y_indmat) * d*d) / sum(y_indmat * d*d)
  return(list(ratio = phi))
}


# local regression on one-on-one paired by p value
pair_by_rank <- function(D, z, y, fun){
  f0_sorted <- get_p(d=D, trt=y, fun=fun)$ratio_all
  fz_sorted <- get_p(mat=z, trt=y, fun=fun)$ratio_all
  N <- length(f0_sorted)
  mat_pair <- matrix(0, nrow=N, ncol=2)
  mat_pair[,1] <- f0_sorted
  mat_pair[,2] <- fz_sorted
  df_pair <- data.frame(data=mat_pair)
  colnames(df_pair) <- c('F0','Fz')
  loess_f <- loess(Fz ~ F0, data=df_pair, span=0.10,
                   control=loess.control(surface="direct"))
  return(list(pair=mat_pair, model=loess_f))
}


# MDS objective term
mds_obj <- function(D, z){
  z_distmat <- as.matrix(dist(z))
  return(sum((D - z_distmat)^2)/2)
}


# confirmatory objective term
conf_obj <- function(y, z, D){
  z_distmat <- as.matrix(dist(z))
  y_indmat <- get_ind_mat(y)
  phi <- sum((1-y_indmat) * D*D) / sum(y_indmat * D*D)
  f_loess <- pair_by_rank(D=D, z=z, y=y, fun=get_phi)$model
  phi_up <- predict(f_loess, phi)  # perform loess for accurate F mapping
  val <- (1-(1+phi_up)*y_indmat) * z_distmat^2
  res <- 0.5 * abs(sum(val))
  return(list(val = res, sign = sign(sum(val))))
}


# gradient term for SGD update, pertaining to a case k=i
grad_obj <- function(N, S, i1, i2, z1, z2, d12, lambda, delta, phi_z, eps){
  grad <- matrix(0, nrow=N, ncol=S)
  grad[i1,] <- (1 - d12/norm(z1-z2, type='2') + lambda*delta*(1-(phi_z+1)*eps)) * (z1 - z2)
  grad[i2,] <- -grad[i1]
  return(grad)
}


# stochastic gradient descent
sgd_cmds <- function(nit = 100, eta = 1e-4, lambda = 0.2, z0, D, y){
  N <- dim(z0)[1]
  S <- dim(z0)[2]
  y_indmat <- get_ind_mat(y)
  phi_o <- sum((1-y_indmat) * D*D) / sum(y_indmat * D*D)
  delta <- conf_obj(y, z0, D)$sign
  z_temp <- z_up <- z0
  F0 <- pseudo_F(d = D, trt = y)$ratio
  
  for(t in 1:nit){
    obj_conf_up <- conf_obj(y, z_up, D)
    obj_mds_up <- mds_obj(D, z_up)
    obj_up <- lambda*obj_conf_up$val + obj_mds_up
    Fz_up <- pseudo_F(mat = z_up, trt = y)$ratio
    print(paste('epoch', t, 
                '  total', sprintf(obj_up, fmt = '%#.3f'), 
                '  mds', sprintf(obj_mds_up, fmt = '%#.3f'), 
                '  conf', sprintf(obj_conf_up$val, fmt = '%#.3f'),
                '  Fz', sprintf(Fz_up, fmt = '%#.2f'),
                '  F0', sprintf(F0, fmt = '%#.2f')))
    if(lambda==0){
      phi_z <- phi_o
    } else {
      f_loess <- pair_by_rank(D=D, z=z_up, y=y, fun=get_phi)$model
      phi_z <- predict(f_loess, phi_o)  # perform loess for accurate F mapping
    }
    inds <- sample(1:N, 2)
    z1 <- z_up[inds[1],]
    z2 <- z_up[inds[2],]
    d12 <- D[inds[1],inds[2]]
    eps12 <- y_indmat[inds[1],inds[2]]
    delta <- conf_obj(y, z_up, D)$sign
    # compute gradient
    grad_obj <- matrix(0, nrow=N, ncol=S)
    grad_obj[inds[1],] <- (1 - d12/norm(z1-z2, type='2') 
                           + lambda*delta*(1-(phi_z+1)*eps12)) * (z1-z2)
    grad_obj[inds[2],] <- -grad_obj[inds[1],]
    z_up <- z_up - eta*grad_obj  # update
  }
  obj_0 <- conf_obj(y, z0, D)$val + lambda*mds_obj(D, z0)
  obj_f <- conf_obj(y, z_up, D)$val + lambda*mds_obj(D, z_up)
  Fz_up <- pseudo_F(mat = z_up, trt = y)$ratio
  return(list(z = z_up, obj_0 = obj_0, obj_f = obj_f, F_z = Fz_up, F_0 = F0))
}


# run
zmds2 <- ordu2$vectors[,1:2]
y2 <- ifelse(site2@sam_data$Treatment == "Pt +", 1, 2)
obj_sgd <- sgd_cmds(nit=100, eta = 1e-3, lambda=20, z0=zmds2, D=distmat2, y=y2)