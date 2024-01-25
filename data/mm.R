# variables dependency on get_dist.R
source("permanova_with_config.R")


# Distance between vector
get_dist_mat <- function(z){
  z_dist <- as.matrix(dist(z))
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
  f_loess <- pair_by_rank(D=D, z=z, y=y, fun=get_phi)$model
  phi_up <- predict(f_loess, phi)  # perform loess for accurate F mapping
  val <- (1-(1+phi_up)*y_indmat) * z_distmat^2
  res <- 0.5 * abs(sum(val))
  return(list(val = res, sign = sign(sum(val))))
}


mm_cmds <- function(nit = 100, conv_crit = 5e-03, lambda = 0.2,
                    z0, D, y){
  N <- dim(z0)[1]
  S <- dim(z0)[2]
  y_indmat <- get_ind_mat(y)
  phi <- sum((1-y_indmat) * D*D) / sum(y_indmat * D*D)
  delta <- conf_obj(y, z0, D)$sign
  z_temp <- z_up <- z0
  # F0 <- pseudo_F(d = D, trt = y)$ratio
  p0 <- get_p(d = D, trt = y)$p
  for(t in 1:nit){
    obj_conf_up <- conf_obj(y, z_up, D)
    obj_mds_up <- mds_obj(D, z_up)
    obj_up <- lambda*obj_conf_up$val + obj_mds_up
    # Fz_up <- pseudo_F(mat = z_up, trt = y)$ratio
    p_up <- get_p(mat = z_up, trt = y)$p
    print(paste('epoch', t, 
                '  total', sprintf(obj_up, fmt = '%#.3f'), 
                '  mds', sprintf(obj_mds_up, fmt = '%#.3f'), 
                '  conf', sprintf(obj_conf_up$val, fmt = '%#.3f'),
                # '  Fz', sprintf(Fz_up, fmt = '%#.2f'),
                # '  F0', sprintf(F0, fmt = '%#.2f')
                '  p_z', sprintf(p_up, fmt = '%#.2f'),
                '  p_0', sprintf(p0, fmt = '%#.2f')
    ))
    
    for(i in 1:N){
      if(lambda==0){
        phi_up <- phi
      } else {
        f_loess <- pair_by_rank(D=D, z=z_up, y=y, fun=get_phi)$model
        phi_up <- predict(f_loess, phi)  # perform loess for accurate F mapping
      }
      
      delta <- conf_obj(y, z_up, D)$sign
      
      z_distmat <- as.matrix(dist(z_up))  # (N,N)
      coeff <- D/z_distmat  # final term in the update
      coeff[is.nan(coeff)] <- 0
      z_diff <- -sweep(x=z_up, MARGIN=2, STATS=as.matrix(z_up[i,]), FUN="-")
      
      z_temp[i,] <- (1+lambda*delta) * (apply(z_up[y!=y[i],], 2, sum)-z_up[i,]) +
        (1-lambda*phi_up*delta) * (apply(z_up[y==y[i],], 2, sum)-z_up[i,]) +
        apply(sweep(x=z_diff, MARGIN=1, STATS=coeff[,i], FUN="*"), 2, sum)
      z_temp[i,] <- z_temp[i,] / (N-1 + 0.5*(N-(N-2)*phi_up)*lambda*delta)
      z_up[i,] <- z_temp[i,]
      if(i %in% seq(from=10, to = N, by=10)){
        print(paste(i, "of", N, "observation updated")) ##to check progress
      }
    }
    # z_up <- z_temp
    # print(z_up)
  }
  
  obj_0 <- conf_obj(y, z0, D)$val + lambda*mds_obj(D, z0)
  obj_f <- conf_obj(y, z_up, D)$val + lambda*mds_obj(D, z_up)
  Fz_up <- pseudo_F(mat = z_up, trt = y)$ratio
  return(list(z = z_up, obj_0 = obj_0, obj_f = obj_f, F_z = Fz_up, F_0 = F0))
}


# run
zmds1 <- ordu1$vectors[,1:2]
zmds2 <- ordu2$vectors[,1:2]
y1 <- ifelse(site1@sam_data$Treatment == "Pt +", 1, 2)
y2 <- ifelse(site2@sam_data$Treatment == "Pt +", 1, 2)
obmm_x <- mm_cmds(nit=15, lambda=0.3, z0=zmds2, D=distmat2, y=y2)  # just an example. replace x with any number you want


# plot
par(mfrow = c(1,2))
ggplot(data.frame(cbind(obmm2$z, (y2s[,1]))), aes(x=Axis.1, y=Axis.2)) + 
  geom_point(aes(color=V3)) + 
  stat_ellipse(level = 0.8, aes(color=V3, group=V3))
plot(obmm0$z, col = y2s[,1])
plot(obmm1.1$z, col = y2s[,1])