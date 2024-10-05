# variables dependency on get_dist.R
source("permanova_with_config.R")


# Distance between vector
get_dist_mat <- function(z){
  z_dist <- as.matrix(dist(z))
  return(z_dist)
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
  return(phi)
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
  # loess_f <- loess(Fz ~ F0, data=df_pair, span=0.10,
  #                  control=loess.control(surface="direct"))
  # return(list(pair=mat_pair, model=loess_f))
  return(list(pair=mat_pair))
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
  list_pair <- pair_by_rank(D=D, z=z, y=y, fun=get_phi)$pair # _0, _z
  ind_phi <- which.min(abs(phi - list_pair[,1]))[1]
  phi_pred <- list_pair[,2][ind_phi]
  # f_loess <- pair_by_rank(D=D, z=z, y=y, fun=get_phi)$model
  # phi_pred <- predict(f_loess, phi)  # perform loess for accurate F mapping
  val <- (1-(1+phi_pred)*y_indmat) * z_distmat^2
  res <- 0.5 * abs(sum(val))
  return(list(val = res, sign = sign(sum(val))))
}


mm_cmds <- function(nit = 100, lambda = 0.2, z0, D, y, dataset = 'example'){
  N <- dim(z0)[1]
  S <- dim(z0)[2]
  y_indmat <- get_ind_mat(y)
  phi <- sum((1-y_indmat) * D*D) / sum(y_indmat * D*D)
  # delta <- conf_obj(y, z0, D)$sign
  z_temp <- z_up <- z0
  # F0 <- pseudo_F(d = D, trt = y)$ratio
  p0 <- get_p(d = D, trt = y)$p
  log_iter_mat <- matrix(0, nrow=0, ncol=6)
  colnames(log_iter_mat) <- c('epoch', 'obj', 'obj_mds', 'obj_confr', 'p_z', 'p_0')
  # obj_prev <- 0
  p_prev <- 1
  
  for(t in 0:nit){
    obj_conf_up <- conf_obj(y, z_up, D)
    obj_mds_up <- mds_obj(D, z_up)
    obj_up <- lambda*obj_conf_up$val + obj_mds_up
    p_up <- get_p(mat = z_up, trt = y)$p
    
    if((abs(p_up-p0) > abs(p_prev-p0)) & (abs(p_prev-p0)<=0.005)){
      print(sprintf('Lambda %.1f ...halt iteration', lambda))
      z_up <- z_prev # revert to prev
      break
    }
    
    print(paste('epoch', t, 
                '  lambda', lambda,
                '  total', sprintf(obj_up, fmt = '%#.2f'), 
                '  mds', sprintf(obj_mds_up, fmt = '%#.2f'), 
                '  conf', sprintf(obj_conf_up$val, fmt = '%#.2f'),
                # '  Fz', sprintf(Fz_up, fmt = '%#.2f'),
                # '  F0', sprintf(F0, fmt = '%#.2f')
                '  p_z', sprintf(p_up, fmt = '%#.3f'),
                '  p_0', sprintf(p0, fmt = '%#.3f')
    ))
    log_iter_mat <- rbind(log_iter_mat, 
                          c(t, obj_up, obj_mds_up, obj_conf_up$val, p_up, p0))
    write.csv(log_iter_mat, sprintf('result/inProg/%s-%.2f-log.csv', dataset, lambda), row.names = FALSE)
    write.csv(z_up, sprintf('result/inProg/%s-%.2f-Z.csv', dataset, lambda), row.names = FALSE)
    
    if(lambda==0){
      phi_pred <- phi
    } else {
      list_pair <- pair_by_rank(D=D, z=z_up, y=y, fun=get_phi)$pair # _0, _z
      ind_phi <- which.min(abs(phi - list_pair[,1]))[1]
      phi_pred <- list_pair[,2][ind_phi]
    }
      
    for(i in 1:N){
      # print(paste('i=', i, Sys.time()))
      # if(lambda==0){
      #   phi_pred <- phi
      # } else {
      #   # f_loess <- pair_by_rank(D=D, z=z_up, y=y, fun=get_phi)$model
      #   # phi_pred <- predict(f_loess, phi)  # perform loess for accurate F mapping
      #   list_pair <- pair_by_rank(D=D, z=z_up, y=y, fun=get_phi)$pair # _0, _z
      #   ind_phi <- which.min(abs(phi - list_pair[,1]))[1]
      #   phi_pred <- list_pair[,2][ind_phi]
      # }
      
      delta <- sign(get_phi(mat=z_up, trt=y) - phi_pred)
      # delta <- conf_obj(y, z_up, D)$sign
      
      z_distmat <- as.matrix(dist(z_up))  # (N,N)
      coeff <- D/z_distmat  # final term in the update
      coeff[is.nan(coeff)] <- 0
      z_diff <- -sweep(x=z_up, MARGIN=2, STATS=as.matrix(z_up[i,]), FUN="-")
      
      z_temp[i,] <- (1+lambda*delta) * (apply(z_up[y!=y[i],], 2, sum)-z_up[i,]) +
        (1-lambda*phi_pred*delta) * (apply(z_up[y==y[i],], 2, sum)-z_up[i,]) +
        apply(sweep(x=z_diff, MARGIN=1, STATS=coeff[,i], FUN="*"), 2, sum)
      z_temp[i,] <- z_temp[i,] / (N-1 + 0.5*(N-(N-2)*phi_pred)*lambda*delta)
      # z_up[i,] <- z_temp[i,]
      # if(i %in% seq(from=10, to = N, by=10)){
      #   print(paste(i, "of", N, "observation updated")) ##to check progress
      # }
    } # end z_temp
    
    z_prev <- z_up
    # obj_prev <- obj_up
    p_prev <- p_up
    z_up <- z_temp
  } # end iteration
  
  obj_0 <- conf_obj(y, z0, D)$val + lambda*mds_obj(D, z0)
  obj_f <- conf_obj(y, z_up, D)$val + lambda*mds_obj(D, z_up)
  Fz_up <- pseudo_F(mat = z_up, trt = y)
  F0 <- pseudo_F(d = D, trt = y)
  write.csv(log_iter_mat, sprintf('result/%s-%.2f-log.csv', dataset, lambda), row.names = FALSE)
  write.csv(z_up, sprintf('result/%s-%.2f-Z.csv', dataset, lambda), row.names=FALSE)
  
  return(list(z = z_up, obj_0 = obj_0, obj_f = obj_f, F_z = Fz_up, F_0 = F0))
}


