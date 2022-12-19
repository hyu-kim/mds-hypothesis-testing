# variables dependency on get_dist.R, data.R
source("permanova_with_config.R")

# Distance between vector
dist_matrix <- function(z){
  N = dim(z)[1]
  z_dist = matrix(0, nrow = N, ncol = N)
  for(i in 1:N){
    for(j in 1:N){
      z_dist[i,j] <- sqrt(sum((z[i,] - z[j,])^2))
    }
  }
  return(z_dist)
}


# Difference between vectors
dist_arr_3d <- function(z){
  N = dim(z)[1]
  S = dim(z)[2]
  z_dist_3d = array(0, dim=c(N,N,S))
  for(i in 1:N){
    for(j in 1:N){
      for(s in 1:S){
        z_dist_3d[i,j,s] <- z[i,s] - z[j,s]
      }
    }
  }
  return(z_dist_3d)
}


# Indicator matrix of labels
ind_matrix <- function(y){
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
  z_dist <- dist_matrix(z)
  return(sum((D - z_dist)^2)/2)
}


# gradient of MDS objective
grad_mds_obj <- function(D, z){
  N <- dim(z)[1]
  S <- dim(z)[2]
  z_dist <- dist_matrix(z)
  z_dist2 <- z_dist + diag(N) # to avoid singularity
  term1 <- 1 - D/z_dist2 # N by N
  term1 <- array(rep(term1, S), dim=c(N,N,S)) # braodcast to size (N,N,S)
  term2 <- dist_arr_3d(z) # (N,N,S)
  res <- 2 * rowSums(aperm(term1*term2, c(1,3,2)), dims=2) # (N, S)
  return(res)
}


# confirmatory objective term with labels
conf_obj <- function(y, z, D){
  N <- length(y)
  z_distmat <- dist_matrix(z)
  y_indmat <- ind_matrix(y)
  phi <- sum((1-y_indmat) * D*D) / sum(y_indmat * D*D)
  res <- (N-2) * sum((1 - (1+phi)*y_indmat) * z_distmat)^2
  return(res)
}


# gradient of confirmatory objective
grad_conf_obj <- function(y, z, D){
  N <- length(y)
  z_distmat <- dist_matrix(z)
  y_indmat <- ind_matrix(y) # (N,N)
  phi <- sum((1-y_indmat) * D*D) / sum(y_indmat * D*D)
  y_ind_3d <- array(rep(y_indmat, S), dim=c(N,N,S))
  z_dist_3d <- dist_arr_3d(z) # (N,N,S)
  term1 <- sum((1 - (1+phi)*y_indmat) * z_distmat)
  term2 <- (1 - (1+phi)*y_ind_3d) * z_dist_3d # (N,N,S)
  term2 <- rowSums(aperm(term2, c(1,3,2)), dims=2) # (N,S)
  res <- 8*(N-2) * term1 * term2
  return(res)
}


# MDS term
gd_mds <- function(nit = 100, eta = 1e-04,
                   z0, distmat){
  N <- dim(z0)[1]
  S <- dim(z0)[2]
  z_cur <- z0 # matrix, N * S
  obj <- mds_obj(D = distmat, z = z_cur)
  for (t in 1:nit){
    print(paste('iteration', t, 
                '  total', sprintf(obj, fmt = '%#.3f')
    ))
    z_cur <- z_cur - eta * grad_mds_obj(D=distmat, z=z_cur)
    obj = mds_obj(D = distmat, z = z_cur)
  }
  obj0 = mds_obj(D = distmat, z = z0)
  
  return(list(z = z_cur, obj0 = obj0, obj_up = obj_up))
}


# MDS + F terms
gd_cmds <- function(nit = 100, eta = 1e-04, conv_crit = 5e-03, lambda = 0.05,
                    z0, distmat, y){
  N <- dim(z0)[1]
  S <- dim(z0)[2]
  z_cur <- z0 #matrix, N * S
  Fz_cur <- pseudo_F(mat = z_cur, trt = y)$pseudoF
  F0 <- pseudo_F(d = distmat, trt = y)$pseudoF
  Finit <- Fz_cur
  obj_total = mds_obj(D = distmat, z = z_cur) + lambda * conf_obj(y=y, z=z_cur, D=distmat)
  for(t in 1:nit){
    obj_total_prev <- obj_total
    print(paste('iteration', t, 
                '  total', sprintf(obj_total, fmt = '%#.3f'), 
                '  Fz', sprintf(Fz_cur, fmt = '%#.3f'),
                '  F0', sprintf(F0, fmt = '%#.3f')
                ))
    #   d_g <- 4 * (N-a)/a * sign(Fz_cur - F0) * (tmp11 * tmp12 - tmp21 * tmp22) * tmp3
    #   z_cur[i,] <- z_cur[i,] - eta * (d_f + lambda * d_g)
    #   Fz_cur <- pseudo_F(mat = z_cur, trt = y)$pseudoF
    # }
    grad_total <- grad_mds_obj(D=distmat, z=z_cur) + lambda * grad_conf_obj(y=y, z=z_cur, D=distmat)
    z_cur <- z_cur - eta * grad_total
    obj_total = mds_obj(D = distmat, z = z_cur) + lambda * conf_obj(y=y, z=z_cur, D=distmat)
    if(abs(obj_total - obj_total_prev) < conv_crit * obj_total){
      print("Converged")
      break
      }
  }
  obj0 = mds_obj(D = distmat, z = z0)
  obj_up = mds_obj(D = distmat, z = z_cur)
  
  return(list(z = z_cur, obj0 = obj0, obj_up = obj_up,
              F0 = F0, Finit = Finit, Fup = Fz_cur))
}


# TEST RUN
y1 <- ifelse(site1@sam_data$Treatment == "Pt +", 1, 2)
y2 <- ifelse(site2@sam_data$Treatment == "Pt +", 1, 2)

tmp <- gd_mds(nit = 100, z0 = zmds1, distmat = as.matrix(dist1))
tmp
plot(tmp$z, col = y2)

# tmp2 <- gd_class(nit = 2, z0 = zmds1, D = as.matrix(dist1), y = y1)
# plot(tmp2$z, col = y1)

tmp3 <- gd_cmds(nit = 30, distmat = as.matrix(dist1), z0 = zmds1, y=y1, lambda = 0.05, eta = 1e-03)
tmp3

tmp4 <- gd_cmds(nit = 15, distmat = as.matrix(dist2), z0 = zmds2, y=y2, lambda = 0.05)
tmp4


par(mfrow = c(1,2))
plot(tmp3$z, col = y1, 
     main = paste("F0=", round(tmp3$F0, 4), ", F_up=", round(tmp3$Fup, 4),
                  "Finit=", round(tmp3$Finit, 4)))
plot(tmp4$z, col = y2,
     main = paste("F0=", round(tmp4$F0, 4), ", F_up=", round(tmp4$Fup, 4),
                  "Finit=", round(tmp4$Finit, 4)))


# EXPORT
wd <- '~/Google Drive/My Drive/Dimension reduction study/2022-11 iss5/results'
sink(paste(wd, "gd_cmds_iter_inf_lamb_5e-2_eta_1e-3.txt", sep='/'))
print(tmp4)
sink()

png(filename=paste(wd, "gd_cmds_iter_inf_lamb_5e-2_eta_1e-3.png", sep='/'))
plot(tmp3$z, col = y1, 
     main = paste("F0=", round(tmp4$F0, 4), 
                  ", F_up=", round(tmp4$Fup, 4),
                  "Finit=", round(tmp4$Finit, 4)))
dev.off()