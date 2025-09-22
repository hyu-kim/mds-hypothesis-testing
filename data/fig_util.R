source('permanova_with_config.R')
## functions
# Distance matrix
get_dist_mat <- function(z){
  z_dist <- dist(z)
  return(z_dist)
}


# Stress-1
get_stress <- function(x_dist, z_dist){
  res <- sqrt(sum((x_dist - z_dist)^2) / sum(z_dist^2))
  return(res)
}


# Raw Stress
get_stress_raw <- function(x_dist, z_dist){
  res <- sum((x_dist - z_dist)^2)
  return(res)
}


# P ratio (do not confuse with earlier version termed p-ratio)
get_p_diff <- function(z, d, y){
  p0 <- get_p(d=d, trt=y)$p
  pemb <- get_p(mat=z, trt=y)$p
  return(abs(p0-pemb)) # ranging [0,1], less variant to data replicates by division
}


# Pearson correlation in distance
get_pearson_corr <- function(x_dist, z_dist){
  res <- cov(x_dist, z_dist) / (sd(x_dist) * sd(z_dist))
  return(res)
}


# p-ratio and correlation in F-ratio
rp_eval2 <- function(nperm=500, dm, y_orig, z_emb, y_emb = NULL){
  n <- length(y_orig)
  if(is.null(y_emb)){
    y_emb <- y_orig
  }
  F0 <- pseudo_F(trt = y_orig, d = dm)
  Femb <- pseudo_F(mat = z_emb, trt = y_emb)
  
  F0_perm <- Femb_perm <- rep(0, nperm)
  for(i in 1:nperm){
    set.seed(i)
    ind_perm <- sample(1:n, size = n, replace = F)
    y0_perm <- y_orig[ind_perm]
    ye_perm <- y_emb[ind_perm]
    
    F0_perm[i] <- pseudo_F(trt = y0_perm, d = dm)
    Femb_perm[i] <- pseudo_F(mat = z_emb, trt = ye_perm)
  }
  
  q = which(Femb == sort(c(Femb_perm, Femb))) /
    which(F0 == sort(c(F0_perm, F0)))
  return(list(F0 = F0, F0_perm = F0_perm,
              Femb = Femb, Femb_perm = Femb_perm,
              rho = cor(F0 - F0_perm, Femb - Femb_perm, method = "pearson"),
              q = q))
}

source('permanova_with_config.R')
## functions
# Distance matrix
get_dist_mat <- function(z){
  z_dist <- dist(z)
  return(z_dist)
}


# Stress-1
get_stress <- function(x_dist, z_dist){
  res <- sqrt(sum((x_dist - z_dist)^2) / sum(z_dist^2))
  return(res)
}


# Raw Stress
get_stress_raw <- function(x_dist, z_dist){
  res <- sum((x_dist - z_dist)^2)
  return(res)
}


# P ratio (do not confuse with earlier version termed p-ratio)
get_p_diff <- function(z, d, y){
  p0 <- get_p(d=d, trt=y)$p
  pemb <- get_p(mat=z, trt=y)$p
  return(abs(p0-pemb)) # ranging [0,1], less variant to data replicates by division
}


# Pearson correlation in distance
get_pearson_corr <- function(x_dist, z_dist){
  res <- cov(x_dist, z_dist) / (sd(x_dist) * sd(z_dist))
  return(res)
}


# p-ratio and correlation in F-ratio
rp_eval2 <- function(nperm=500, dm, y_orig, z_emb, y_emb = NULL){
  n <- length(y_orig)
  if(is.null(y_emb)){
    y_emb <- y_orig
  }
  F0 <- pseudo_F(trt = y_orig, d = dm)
  Femb <- pseudo_F(mat = z_emb, trt = y_emb)
  
  F0_perm <- Femb_perm <- rep(0, nperm)
  for(i in 1:nperm){
    set.seed(i)
    ind_perm <- sample(1:n, size = n, replace = F)
    y0_perm <- y_orig[ind_perm]
    ye_perm <- y_emb[ind_perm]
    
    F0_perm[i] <- pseudo_F(trt = y0_perm, d = dm)
    Femb_perm[i] <- pseudo_F(mat = z_emb, trt = ye_perm)
  }
  
  q = which(Femb == sort(c(Femb_perm, Femb))) /
    which(F0 == sort(c(F0_perm, F0)))
  return(list(F0 = F0, F0_perm = F0_perm,
              Femb = Femb, Femb_perm = Femb_perm,
              rho = cor(F0 - F0_perm, Femb - Femb_perm, method = "pearson"),
              q = q))
}