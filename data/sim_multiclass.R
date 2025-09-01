## Simulation for trinary dataset
# import package
devtools::install_github("biobakery/SparseDOSSA2")
library(vegan)
library(SparseDOSSA2)
library(MASS)
library(parallel)
source('mm.R')


## Load stool, parameterize, perturb, and export to a template
get_template <- function(data_str = "Stool"){
  
  # obtain feature params from stool dataset
  stool <- SparseDOSSA2(template = data_str, n_sample = 100, n_feature = 332, verbose = TRUE, new_features = FALSE)
  params <- stool$params$feature_param
  params <- cbind(params, mean=0)
  params <- cbind(params, var=0)
  params[,'mean'] <- exp(params[,'mu'] + params[,'sigma']^2/2)
  params[,'var'] <- (exp(params[,'sigma']^2)-1) * exp(2*params[,'mu'] + params[,'sigma']^2)
  params <- cbind(params, sigma_adj=0)
  params[,'sigma_adj'] <- sqrt(
    (1-params[,'pi0']) * params[,'var'] + params[,'pi0']*(1-params[,'pi0']) * params[,'mean']^2
  )
  
  # select features for perturbation
  params_df <- as.data.frame(params)
  n_ind_diff <- 5
  
  # # Two features with high variance and for zero mean difference
  # ind_same <- order(params_df$sigma_adj, decreasing = TRUE)[1:2]
  # features1_list <- rownames(params_df)[ind_same]
  # n features with low variance and for high mean difference
  ind_diff <- order(params_df$sigma_adj, decreasing = FALSE)[1:n_ind_diff]
  features2_list <- rownames(params_df)[ind_diff]
  
  # perturb features
  # params_df$sigma[ind_diff] <- params_df$sigma[ind_diff] * 0.5
  params_df3 <- params_df2 <- params_df1 <- params_df
  
  i <- seq_len(n_ind_diff)
  offset <- (5 - (i - 1) %/% (n_ind_diff %/% 5)) * (-1)^(i+1)
  
  params_df2$mu[ind_diff] <- params_df$mu[ind_diff] + offset
  params_df3$mu[ind_diff] <- params_df$mu[ind_diff] + 2*offset
  

  # templates based on params
  template3 <- template2 <- template1 <- stool$template
  
  template1$EM_fit$fit$mu[] <- params_df1$mu
  template2$EM_fit$fit$mu[] <- params_df2$mu
  template3$EM_fit$fit$mu[] <- params_df3$mu
  
  template1$EM_fit$fit$sigma[] <- params_df1$sigma
  template2$EM_fit$fit$sigma[] <- params_df2$sigma
  template3$EM_fit$fit$sigma[] <- params_df3$sigma
  
  templates <- list()
  templates$"1" <- template1
  templates$"2" <- template2
  templates$"3" <- template3
  
  return(templates)
}


## simulate new dataset
get_sim_data_mat <- function(templates, N = 25){
  stool1 <- SparseDOSSA2(template = templates$"1", n_sample = N, n_feature = 332, 
                         verbose = FALSE, new_features = FALSE)
  stool2 <- SparseDOSSA2(template = templates$"2", n_sample = N, n_feature = 332, 
                         verbose = FALSE, new_features = FALSE)
  stool3 <- SparseDOSSA2(template = templates$"3", n_sample = N, n_feature = 332, 
                         verbose = FALSE, new_features = FALSE)
  
  mat1 <- stool1$simulated_matrices$rel
  mat2 <- stool2$simulated_matrices$rel
  mat3 <- stool3$simulated_matrices$rel
  
  colnames(mat1) <- paste(colnames(mat1), "-1", sep="")
  colnames(mat2) <- paste(colnames(mat2), "-2", sep="")
  colnames(mat3) <- paste(colnames(mat3), "-3", sep="")
  
  mat <- cbind(mat1, mat2, mat3)
  
  return(mat)
}


## Perform statistical tests and print results
get_p_values <- function(N, mat){
  group <- factor(rep(c(0, 1, 2), each = N))
  
  # PERMANOVA, original
  dist <- vegdist(t(mat), method="bray")
  result <- adonis2(dist ~ group)
  p0 <- result$`Pr(>F)`[1]
  
  # PERMANOVA, 2D
  pcoa <- cmdscale(dist, eig = TRUE)
  dist2 <- dist(pcoa$points[,1:2], method='euclidean')
  result2 <- adonis2(dist2 ~ group)
  pz <- result2$`Pr(>F)`[1]
  
  print(sprintf("Original P-value = %g, 2-D P-value = %g", p0, pz))
  
  return(list(p0 = p0, pz = pz))
}


## Run until P value does not match
N <- 25 # set half of total size
templates <- get_template()

p0 <- 1
pz <- 1

while (p0>0.05 | pz<0.5) {
  print(sprintf('Data of size %d does not outlie. Generating new..', 3*N))
  mat <- get_sim_data_mat(templates, N = N)
  p_list <- get_p_values(N, mat)
  p0 <- p_list$p0
  pz <- p_list$pz
}

print('done.')


## visualize PCoA
dist <- vegdist(t(mat), method="bray")
pcoa <- cmdscale(dist, eig = TRUE)
graphics.off()
plot(pcoa$points[,1], pcoa$points[,2], col='white')
points(pcoa$points[1:N,1], pcoa$points[1:N,2], col='red')
points(pcoa$points[(N+1):(2*N),1], pcoa$points[(N+1):(2*N),2], col='blue')
points(pcoa$points[(2*N+1):(3*N),1], pcoa$points[(2*N+1):(3*N),2], col='black')


## save data
y <- rep(c(0, 1, 2), each = N)
saveRDS(mat, "result/Multiclass/sim4d_rev-data.rds")
saveRDS(dist, "result/Multiclass/sim4d_rev-dist.rds")
write.csv(y, 'result/Multiclass/sim4d_dist-Y.csv', row.names=FALSE)


## run FMDS
res <- mm_cmds(nit=50, lambda=0, z0=pcoa$points, D=as.matrix(dist), y=y, dataset = 'sim4d_rev')
res <- mm_cmds(nit=100, lambda=0.5, z0=pcoa$points, D=as.matrix(dist), y=y, dataset = 'sim4d_rev')
res <- mm_cmds(nit=100, lambda=1, z0=pcoa$points, D=as.matrix(dist), y=y, dataset = 'sim4d_rev')