# variables dependency on get_dist.R, data.R
source("permanova_with_config.R")

# MDS objective
mds_obj <- function(D, z){
  N = dim(D)[1]
  S = dim(z)[2]
  z_dist = matrix(0, nrow = N, ncol = N)
  for(i in 1:N){
    for(j in 1:N){
      z_dist[i,j] <- sqrt(sum((z[i,] - z[j,])^2))
    }
  }
  return(sum((D - z_dist)^2)/2)
}


# -----
# MDS term
gd_mds <- function(nit = 1000, eta = 1e-04, conv_lim = 1e-05,
                   D, S = 2, z0){
  # fixed step size (eta)
  N = dim(D)[1]
  z_cur <- z0 #matrix, N * S
  for(t in 1:nit){
    print(paste('iteration',t))
    for(i in 1:N){
      # f
      d_f <- rep(0, S)
      for(j in setdiff(1:N, i)){
        d_f <- d_f + 2 * (1 - D[i,j]/sqrt(sum((z_cur[i,] - z_cur[j,])^2))) * 
          (z_cur[i,] - z_cur[j,])
      }
      z_cur[i,] <- z_cur[i,] - eta * d_f
    }
    z_up <- z_cur
    # if(sqrt(sum((z_up - z_cur)^2)) < conv_lim){
    #   print("Converged")
    #   break
    #### convergence at every i?
    # }
  }
  
  obj0 = mds_obj(D = D, z = z0)
  obj_up = mds_obj(D = D, z = z_up)
  return(list(z = z_up, obj0 = obj0, obj_up = obj_up))
}


# -----
# F statistic term, L1 norm
gd_class <- function(nit = 1000, eta = 1e-04, z0, conv_lim = 1e-05,
                     D, S = 2, y){
  # fixed step size (eta)
  N = nrow(z0)
  a <- length(table(y))
  
  z_cur <- z0 #matrix, N * S
  Fz_cur <- pseudo_F(mat = z_cur, trt = y)$pseudoF
  F0 <- pseudo_F(d = D, trt = y)$pseudoF
  Finit <- Fz_cur
  for(t in 1:nit){
    for(i in 1:N){
      # f
      # d_g <- rep(0, S)
      tmp11 <- tmp22 <- rep(0, S)
      tmp12 <- tmp21 <- 0
      for(j in 1:N){
        tmp11 <- tmp11 + z_cur[i,] - z_cur[j,]
        tmp22 <- tmp22 + ifelse(y[i] == y[j], z_cur[i,] - z_cur[j,], 0)
        for(i_ in 1:N){
          tmp21 <- tmp21 + sum((z_cur[i_,] - z_cur[j,])^2)
          tmp12 <- tmp12 + ifelse(y[i_] == y[j], 
                                  sum((z_cur[i_,] - z_cur[j,])^2), 0)
        }
      }
      d_g <- 4 * (N-a)/a * sign(Fz_cur - F0) * (tmp11 * tmp12 - tmp21 * tmp22)
      z_cur[i,] <- z_cur[i,] - eta * d_g
      Fz_cur <- pseudo_F(mat = z_cur, trt = y)$pseudoF
    }
    z_up <- z_cur
    # if(sqrt(sum((z_up - z_cur)^2)) < conv_lim){
    #   print("Converged")
    #   break
    #### convergence at every i?
    # }
  }
  return(list(z = z_up, obj0 = abs(F0 - Finit), obj_up = abs(F0 - Fz_cur),
              F0 = F0, Finit = Finit, Fup = Fz_cur))
}


# -----
# MDS + F terms
gd_cmds <- function(nit = 1000, eta = 1e-04, conv_crit = 5e-03, lambda = 0.05,
                    z0, D, S = 2, y){
  # fixed step size (eta)
  N = dim(D)[1]
  a <- length(table(y))
  z_cur <- z0 #matrix, N * S
  Fz_cur <- pseudo_F(mat = z_cur, trt = y)$pseudoF
  F0 <- pseudo_F(d = D, trt = y)$pseudoF
  Finit <- Fz_cur
  obj_total = mds_obj(D = D, z = z_cur) + lambda*abs(F0 - Fz_cur)
  for(t in 1:nit){
    obj_total_prev <- obj_total
    print(paste('iteration', t, 
                '  total', sprintf(obj_total, fmt = '%#.3f'), 
                '  Fz', sprintf(Fz_cur, fmt = '%#.3f')))
    for(i in 1:N){
      d_f <- rep(0, S)
      tmp11 <- tmp22 <- rep(0, S)
      tmp12 <- tmp21 <- 0
      for(j in 1:N){
        if(j != i) {
          d_f <- d_f + 2 * (1 - D[i,j]/sqrt(sum((z_cur[i,] - z_cur[j,])^2))) * 
            (z_cur[i,] - z_cur[j,])
        }
        tmp11 <- tmp11 + z_cur[i,] - z_cur[j,]
        tmp22 <- tmp22 + as.numeric(ifelse(y[i] == y[j], z_cur[i,] - z_cur[j,], rep(0, S)))
        for(i_ in 1:N){
          tmp21 <- tmp21 + sum((z_cur[i_,] - z_cur[j,])^2)
          tmp12 <- tmp12 + ifelse(y[i_] == y[j], 
                                  sum((z_cur[i_,] - z_cur[j,])^2), 0)
        }
      }
      tmp3 <- tmp12^(-2)
      d_g <- 4 * (N-a)/a * sign(Fz_cur - F0) * (tmp11 * tmp12 - tmp21 * tmp22) * tmp3
      z_cur[i,] <- z_cur[i,] - eta * (d_f + lambda * d_g)
      Fz_cur <- pseudo_F(mat = z_cur, trt = y)$pseudoF
    }
    obj_total = mds_obj(D = D, z = z_cur) + lambda*abs(F0 - Fz_cur)
    if(abs(obj_total - obj_total_prev) < conv_crit * obj_total){
      print("Converged")
      break
      }
  }
  obj0 = mds_obj(D = D, z = z0)
  obj_up = mds_obj(D = D, z = z_cur)
  
  return(list(z = z_cur, obj0 = obj0, obj_up = obj_up,
              F0 = F0, Finit = Finit, Fup = Fz_cur))
}


# TEST RUN
y1 <- ifelse(site1@sam_data$Treatment == "Pt +", 1, 2)
y2 <- ifelse(site2@sam_data$Treatment == "Pt +", 1, 2)

tmp <- gd_mds(nit = 100, D = as.matrix(dist1), z0 = zmds1)
tmp
plot(tmp$z, col = y2)

tmp2 <- gd_class(nit = 2, z0 = zmds1, D = as.matrix(dist1), y = y1)
plot(tmp2$z, col = y1)

tmp3 <- gd_cmds(nit = 30, D = as.matrix(dist1), z0 = zmds1, y=y1, lambda = 0.05, eta = 1e-03)
tmp3

tmp4 <- gd_cmds(nit = 15, D = as.matrix(dist2), z0 = zmds2, y=y2, lambda = 0.05)
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