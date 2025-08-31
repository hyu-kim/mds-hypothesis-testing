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
  stool <- SparseDOSSA2(template = data_str, n_sample = 50, n_feature = 332, verbose = TRUE, new_features = FALSE)
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
  
  # Two features with high variance and for zero mean difference
  ind1 <- order(params_df$sigma_adj, decreasing = TRUE)[1:2]
  features1_list <- rownames(params_df)[ind1]
  # One feature with low variance and for high mean difference
  ind2 <- order(params_df$sigma_adj, decreasing = FALSE)[1]
  features2_list <- rownames(params_df)[ind2]
  
  # perturb features
  params_df2 <- params_df1 <- params_df
  params_df1$mu[ind2] <- params_df1$mu[ind2] + 5
  params_df2$mu[ind2] <- params_df2$mu[ind2] - 5
  params_df1$mu[ind1[2]] <- params_df1$mu[ind1[2]] - 0.2
  params_df2$mu[ind1[2]] <- params_df2$mu[ind1[2]] + 0.2
  
  # templates based on params
  template2 <- template1 <- stool$template
  template1$EM_fit$fit$mu[] <- params_df1$mu
  template2$EM_fit$fit$mu[] <- params_df2$mu
  template1$EM_fit$fit$sigma[] <- params_df1$sigma
  template2$EM_fit$fit$sigma[] <- params_df2$sigma
  
  templates <- list()
  templates$"1" <- template1
  templates$"2" <- template2
  
  return(templates)
}


set.seed(150)
sim_data <- list(
  data = rbind(
    mvrnorm(50, c(0,0,0,0), 
            Sigma = matrix(c(5,0,0,0, 0,5,0,0, 0,0,1,0, 0,0,0,1), nrow=4)),
    mvrnorm(50, c(0,0,2,0), 
            Sigma = matrix(c(5,0,0,0, 0,5,0,0, 0,0,1,0, 0,0,0,1), nrow=4)),
    mvrnorm(50, c(0,0,1,sqrt(3)), 
            Sigma = matrix(c(5,0,0,0, 0,5,0,0, 0,0,1,0, 0,0,0,1), nrow=4))
    )
)
sim_data$data <- data.frame(sim_data$data)
sim_data$Y <- c(rep(1, 50), rep(2, 50), rep(3, 50))

sim_data$distmat <- as.matrix(dist(sim_data$data[, 1:4]))

z0 <- cmdscale(dist(sim_data$data[,1:4]), k = 2)

write.csv(sim_data$data, 'result/Multiclass/sim4d_2-data.csv', row.names=FALSE)
write.csv(sim_data$Y, 'result/Multiclass/sim4d_2-Y.csv', row.names=FALSE)
write.csv(z0, 'result/Multiclass/sim4d_2-mds-Z.csv', row.names=FALSE)

res <- mm_cmds(nit=100, lambda=0.5, z0=z0, D=sim_data$distmat, y=sim_data$Y, dataset = 'sim4d_2')

