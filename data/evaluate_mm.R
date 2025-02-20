source('mm.R')

mm_cmds_dev <- function(nit = 50, lambda = 0.4, z0, D, y, dataset = 'example'){
  N <- dim(z0)[1]
  S <- dim(z0)[2]
  a <- length(unique(y))
  y_indmat <- get_ind_mat(y)
  f_ratio <- pseudo_F(d = D, trt = y)
  z_temp <- z_up <- z0
  p0 <- get_p(d = D, trt = y)$p
  log_iter_mat <- matrix(0, nrow=0, ncol=8)
  colnames(log_iter_mat) <- c('epoch', 'F_z', 'F_pred', 'obj', 'obj_mds', 'obj_confr', 'p_z', 'p_0')
  # obj_prev <- 0
  p_prev <- 1
  
  for(t in 0:nit){
    obj_conf_up <- conf_obj(y, z_up, D)
    obj_mds_up <- mds_obj(D, z_up)
    obj_up <- lambda*obj_conf_up$val + obj_mds_up
    p_up <- get_p(mat = z_up, trt = y)$p
    
    if((abs(p_up-p0) > abs(p_prev-p0)) & (abs(p_prev-p0)<=0.05)){
      print(sprintf('Lambda %.2f ...halt iteration', lambda))
      z_up <- z_prev # revert to prev
      break
    }
    
    if(lambda==0){
      f_ratio_pred <- f_ratio
    } else {
      list_pair <- pair_by_rank(D=D, z=z_up, y=y, fun=pseudo_F)$pair # _0, _z
      ind_f_ratio <- which.min(abs(f_ratio - list_pair[,1]))[1]
      f_ratio_pred <- list_pair[,2][ind_f_ratio]
    }
    
    f_ratio_z <- pseudo_F(mat=z_up, trt=y)
    delta <- sign(f_ratio_z - f_ratio_pred)
    z_distmat <- as.matrix(dist(z_up))
    f_diff_denominator <- 2 * sum(y_indmat * z_distmat^2)
    f_diff_nominator <- sum((1 - 2*y_indmat*(1+f_ratio_pred/(N-2))) * z_distmat^2)
    
    print(paste('epoch', t, 
                '  lambda', lambda,
                '  Fz', sprintf(f_ratio_z, fmt = '%#.2f'),
                '  F_pred', sprintf(f_ratio_pred, fmt = '%#.2f'),
                '  denom', sprintf(f_diff_denominator, fmt = '%#.2f'),
                '  nom', sprintf(f_diff_nominator, fmt = '%#.2f'),
                # '  total', sprintf(obj_up, fmt = '%#.2f'), 
                '  mds', sprintf(obj_mds_up, fmt = '%#.2f'), 
                '  conf', sprintf(obj_conf_up$val, fmt = '%#.2f'),
                '  p_z', sprintf(p_up, fmt = '%#.3f'),
                '  p_0', sprintf(p0, fmt = '%#.3f')
    ))
    log_iter_mat <- rbind(log_iter_mat, 
                          c(t, 
                            f_ratio_z, f_ratio_pred, 
                            f_diff_denominator,
                            # obj_up, 
                            obj_mds_up, obj_conf_up$val, 
                            p_up, p0))
    # write.csv(log_iter_mat, sprintf('result/CompareTerms/%s-fmds-%.2f-%02d-log.csv', dataset, lambda, t), row.names = FALSE)
    write.csv(z_up, sprintf('result/CompareTerms/%s-fmds-%.2f-%02d-Z.csv', dataset, lambda, t), row.names = FALSE)
    
    for(i in 1:N){
      z_distmat <- as.matrix(dist(z_up))  # (N,N)
      coeff <- D/z_distmat  # final term in the update
      coeff[is.nan(coeff)] <- 0
      z_diff <- -sweep(x=z_up, MARGIN=2, STATS=as.matrix(z_up[i,]), FUN="-")
      
      z_temp[i,] <- (1+lambda*delta) * (apply(z_up[y!=y[i],], 2, sum)) +
        (1-lambda*delta*(1+2*f_ratio_pred/(N-2))) * (apply(z_up[y==y[i],], 2, sum)) +
        apply(sweep(x=z_diff, MARGIN=1, STATS=coeff[,i], FUN="*"), 2, sum)
      
      z_temp[i,] <- z_temp[i,] / (N - N*lambda*delta*f_ratio_pred/(N-2))
    } # end z_temp
    
    z_prev <- z_up
    p_prev <- p_up
    z_up <- z_temp
  } # end iteration
  
  obj_0 <- conf_obj(y, z0, D)$val + lambda*mds_obj(D, z0)
  obj_f <- conf_obj(y, z_up, D)$val + lambda*mds_obj(D, z_up)
  Fz_up <- pseudo_F(mat = z_up, trt = y)
  F0 <- pseudo_F(d = D, trt = y)
  write.csv(log_iter_mat, sprintf('result/CompareTerms/%s-fmds-%.2f-final-log.csv', 
                                  dataset, lambda), row.names = FALSE)
  write.csv(z_up, sprintf('result/CompareTerms/%s-fmds-%.2f-final-Z.csv', 
                          dataset, lambda), row.names=FALSE)
  
  return(list(z = z_up, obj_0 = obj_0, obj_f = obj_f, F_z = Fz_up, F_0 = F0))
}


# test example
sim_data <- as.matrix(read.csv('result/Evaluation/sim_1-data.csv'))
y_sim <- as.matrix(read.csv('result/Evaluation/sim_1-Y.csv'))
z0 <- read.csv('result/Evaluation/sim_1-mds-Z.csv')
lambda = 0.4

res <- mm_cmds_dev(nit=50, lambda=lambda, z0=z0, D=as.matrix(dist(sim_data)), y = y_sim, dataset = 'sim')