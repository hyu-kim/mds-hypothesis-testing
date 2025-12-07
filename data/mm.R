# variables dependency on get_dist.R
source("permanova_with_config.R")


# Distance between vector
get_dist_mat <- function(z){
  z_dist <- as.matrix(dist(z))
  return(z_dist)
}


# local regression on one-on-one paired by p value
pair_by_rank <- function(D, z, y, fun){
  f_sorted <- get_f_pair(z=z, d=D, trt=y)
  return(f_sorted)
}


# get lambda from models
get_lambda <- function(p=p_up, p_diff=p_up-p0, model_lambda=model_out){
  lambda <- predict(model_lambda, newdata=data.frame(p_z=p, p_z_diff=p_diff))
  if(is.na(lambda)){
    if (p_diff>0.11) lambda <- 0.5
    else lambda <- 0.08
  }
  return(lambda)
}


# Raw Stress
mds_obj <- function(D, z){
  N <- nrow(z)
  z_distmat <- get_dist_mat(z)
  return(sum((D - z_distmat)^2)/N) # normalize - from scaling analysis
}


# MDS objective normalized (stress-1)
# mds_obj_norm <- function(D, z){
#   z_distmat <- get_dist_mat(z)
#   return(sum((D - z_distmat)^2)/sum(z_distmat^2))
# }


# confirmatory objective term with labels
conf_obj <- function(y, z, D){
  N <- length(y)
  a <- length(unique(y))
  z_distmat <- get_dist_mat(z)
  y_indmat <- get_ind_mat(y)
  f_ratio <- pseudo_F(d = D, trt = y)
  list_pair <- pair_by_rank(D=D, z=z, y=y, fun=pseudo_F) # _0, _z
  ind_f_ratio <- which.min(abs(f_ratio - list_pair[,1]))[1]
  f_ratio_pred <- list_pair[,2][ind_f_ratio]
  val <- (1 - a * y_indmat * (1 + f_ratio_pred*(a-1)/(N-a))) * z_distmat^2
  res <- abs(sum(val))
  return(list(val = res, sign = sign(sum(val))))
}


mm_cmds <- function(nit = 100, lambda = 0.2, z0, D, y, dataset = 'example',
                    threshold = 0.1, folder = "result/Revision2_Dec2025"){
  if(is.null(D)){
    D <- getDistMat(X)
  } else {
    D <- as.matrix(D)
  }
  if(is.null(z0)){
    z0 <- cmdscale(d = D)
  }
  N <- dim(z0)[1]
  S <- dim(z0)[2]
  a <- length(unique(y))
  y_indmat <- get_ind_mat(y)
  f_ratio <- pseudo_F(mat=z, d = D, trt = y)
  z_temp <- z_up <- z0
  p0 <- get_p(d = D, trt = y)$p
  log_iter_mat <- matrix(0, nrow=0, ncol=6)
  colnames(log_iter_mat) <- c('epoch', 'obj', 'obj_mds', 'obj_confr', 'p_z', 'p_0')
  # obj_prev <- 0
  p_prev <- 1
  
  for(t in 0:nit){
    p_up <- get_p(mat = z_up, trt = y)$p
    
    if((abs(p_up-p0) > abs(p_prev-p0)) & (abs(p_prev-p0)<=threshold)){
      print(sprintf('Lambda %.2f ...halt iteration', lambda))
      z_up <- z_prev # revert to prev
      break
    }
    
    # if(lambda==0){
    #   f_ratio_pred <- f_ratio
    # } else {
    list_pair <- pair_by_rank(D=D, z=z_up, y=y, fun=pseudo_F) # _0, _z
    ind_f_ratio <- which.min(abs(f_ratio - list_pair[,1]))[1]
    f_ratio_pred <- list_pair[,2][ind_f_ratio]
    # }
      
    z_distmat <- as.matrix(dist(z_up))
    f_diff_nominator <- sum((1 - a * y_indmat * (1+f_ratio_pred*(a-1)/(N-a))) * z_distmat^2)
    delta <- sign(f_diff_nominator)
    obj_conf <- abs(f_diff_nominator)
    obj_mds <- mds_obj(D, z_up)
    # obj_mds_norm <- mds_obj_norm(D, z_up)
    obj <- lambda*obj_conf + obj_mds
    
    message(paste('epoch', t, 
                '  lambda', lambda,
                # '  str-1', sprintf(obj_mds_norm, fmt = '%#.2f'), 
                '  total', sprintf(obj, fmt = '%#.2f'), 
                '  mds', sprintf(obj_mds, fmt = '%#.2f'), 
                '  conf', sprintf(obj_conf, fmt = '%#.2f'),
                '  p_z', sprintf(p_up, fmt = '%#.3f'),
                '  p_0', sprintf(p0, fmt = '%#.3f')
    ))
    log_iter_mat <- rbind(log_iter_mat, 
                          c(t, obj, obj_mds, obj_conf, p_up, p0))
    
    # update rule
    for(i in 1:N){
      z_distmat <- as.matrix(dist(z_up))  # (N,N)
      coeff <- D/z_distmat  # final term in the update
      coeff[is.nan(coeff)] <- 0
      z_diff <- -sweep(x=z_up, MARGIN=2, STATS=as.matrix(z_up[i,]), FUN="-")
      
      z_temp[i,] <- (1 + lambda*delta) * (apply(z_up[y!=y[i],], 2, sum)) +
        (1 - lambda*delta*(a-1 + f_ratio_pred*a*(a-1)/(N-a))) * (apply(z_up[y==y[i],], 2, sum)) +
        1 * apply(sweep(x=z_diff, MARGIN=1, STATS=coeff[,i], FUN="*"), 2, sum)
      
      z_temp[i,] <- z_temp[i,] / (N - N * lambda*delta*f_ratio_pred/(N-a))
    } # end z_temp
    
    z_prev <- z_up
    # obj_prev <- obj_up
    p_prev <- p_up
    z_up <- z_temp
  } # end iteration
  
  Fz_up <- pseudo_F(mat = z_up, trt = y)
  F0 <- pseudo_F(d = D, trt = y)
  write.csv(log_iter_mat, sprintf('%s/%s-fmds-%.2f-log.csv', 
                                  folder, dataset, lambda), row.names = FALSE)
  write.csv(z_up, sprintf('%s/%s-fmds-%.2f-Z.csv', 
                          folder, dataset, lambda), row.names=FALSE)
  
  return(list(z = z_up, 
              # obj_0 = obj_0, obj_f = obj_f, 
              F_z = Fz_up, F_0 = F0))
}