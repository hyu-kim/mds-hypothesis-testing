# variables dependency on get_dist.R
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
  F0 <- pseudo_F(d = D, trt = y)$pseudoF
  for(t in 1:nit){
    obj_conf_up <- conf_obj(y, z_up, D)
    obj_mds_up <- mds_obj(D, z_up)
    obj_up <- lambda*obj_conf_up$val + obj_mds_up
    Fz_up <- pseudo_F(mat = z_up, trt = y)$pseudoF
    print(paste('epoch', t, 
                '  total', sprintf(obj_up, fmt = '%#.3f'), 
                '  mds', sprintf(obj_mds_up, fmt = '%#.3f'), 
                '  conf', sprintf(obj_conf_up$val, fmt = '%#.3f'),
                '  Fz', sprintf(Fz_up, fmt = '%#.2f'),
                '  F0', sprintf(F0, fmt = '%#.2f')
    ))
    
    delta <- obj_conf_up$sign  # N
    z_distmat <- get_dist_mat(z_up)  # (N,N)
    coeff <- D/z_distmat  # final term in the update
    coeff[is.nan(coeff)] <- 0
    for(i in 1:N){
      z_diff <- -sweep(x=z_up, MARGIN=2, STATS=as.matrix(z_up[i,]), FUN="-")
      z_temp[i,] <- (1+lambda*delta[i]) * (apply(z_up[y!=y[i],], 2, sum)-z_up[i,]) +
        (1-lambda*phi*delta[i]) * (apply(z_up[y==y[i],], 2, sum)-z_up[i,]) +
        apply(sweep(x=z_diff, MARGIN=1, STATS=coeff[,i], FUN="*"), 2, sum)
      z_temp[i,] <- z_temp[i,] / (N-1 + 0.5*(N-(N-2)*phi)*lambda*delta[i])
      z_up[i,] <- z_temp[i,]
    }
    # z_up <- z_temp
    # print(z_up)
  }
  
  obj_0 <- conf_obj(y, z0, D)$val + lambda*mds_obj(D, z0)
  obj_f <- conf_obj(y, z_up, D)$val + lambda*mds_obj(D, z_up)
  Fz_up <- pseudo_F(mat = z_up, trt = y)$pseudoF
  return(list(z = z_up, obj_0 = obj_0, obj_f = obj_f, F_z = Fz_up, F_0 = F0))
}


# run
zmds1 <- ordu1$vectors[,1:2]
zmds2 <- ordu2$vectors[,1:2]
y1 <- ifelse(site1@sam_data$Treatment == "Pt +", 1, 2)
y2 <- ifelse(site2@sam_data$Treatment == "Pt +", 1, 2)
obmm3 <- mm_cmds(nit=200, lambda=0.3, z0=zmds2, D=distmat2, y=y2s[,1])


# plot
par(mfrow = c(1,2))
plot(obmm0$z, col = y2s[,1])
plot(obmm3$z, col = y2s[,1])