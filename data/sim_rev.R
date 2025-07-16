# import package
devtools::install_github("biobakery/SparseDOSSA2")
library(vegan)
library(SparseDOSSA2)

# obtain feature params from stool dataset
stool <- SparseDOSSA2(template = "Stool", n_sample = 50, n_feature = 332, verbose = TRUE, new_features = FALSE)
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

# Four features with high variance and for zero mean difference
ind1 <- order(params_df$sigma_adj, decreasing = TRUE)[1:2]
features1_list <- rownames(params_df)[ind1]
# Two features with low variance and for high mean difference
ind2 <- order(params_df$sigma_adj, decreasing = FALSE)[1]
features2_list <- rownames(params_df)[ind2]
# Rest features for low mean difference
ind3 <- setdiff(1:nrow(params_df), c(ind1, ind2))
features3_list <- setdiff(rownames(params_df), c(features1_list, features2_list))

# perturb features
params_df2 <- params_df1 <- params_df
params_df1$mu[ind2] <- params_df1$mu[ind2] + 5
params_df2$mu[ind2] <- params_df2$mu[ind2] - 5
# params_df1$mu[ind3] <- params_df1$mu[ind3] - 0.1
# params_df2$mu[ind3] <- params_df2$mu[ind3] + 0.1
params_df1$mu[ind1[2]] <- params_df1$mu[ind1[2]] - 0.2
params_df2$mu[ind1[2]] <- params_df2$mu[ind1[2]] + 0.2
# params_df1$sigma[ind1] <- params_df1$sigma[ind1] * 1.5
# params_df2$sigma[ind1] <- params_df2$sigma[ind1] * 1.5
# params_df1$sigma[ind2] <- params_df1$sigma[ind2] * 0.5
# params_df2$sigma[ind2] <- params_df2$sigma[ind2] * 0.5

# templates based on params
template2 <- template1 <- stool$template
template1$EM_fit$fit$mu[] <- params_df1$mu
template2$EM_fit$fit$mu[] <- params_df2$mu
template1$EM_fit$fit$sigma[] <- params_df1$sigma
template2$EM_fit$fit$sigma[] <- params_df2$sigma

# load saved template
template <- readRDS("result/HyperparameterStudy/sim_rev-template.Rds")
template1 <- template$"1"
template2 <- template$"2"

# simulate new dataset
N <- 500
stool1 <- SparseDOSSA2(template = template1, n_sample = N, n_feature = 332, verbose = TRUE, new_features = FALSE)
stool2 <- SparseDOSSA2(template = template2, n_sample = N, n_feature = 332, verbose = TRUE, new_features = FALSE)
mat1 <- stool1$simulated_matrices$rel
mat2 <- stool2$simulated_matrices$rel
colnames(mat2) <- paste(colnames(mat2), "-2", sep="")
mat <- cbind(mat1, mat2)


# load saved result if necessary
mat <- readRDS(sprintf("result/ScalingStudy/sim_rev_1-N%d-data.Rds", 2*N))
group <- read.csv(sprintf("result/ScalingStudy/sim_rev-N%d-Y.csv", 2*N))

# visualize PCoA
dist <- vegdist(t(mat), method="bray")
# dist <- as.matrix(dist)
pcoa <- cmdscale(dist, eig = TRUE)
graphics.off()
plot(pcoa$points[,1], pcoa$points[,2], col='white')
points(pcoa$points[1:N,1], pcoa$points[1:N,2], col='red')
points(pcoa$points[(N+1):(2*N),1], pcoa$points[(N+1):(2*N),2], col='blue')

# PERMANOVA, original
group <- factor(rep(c(0, 1), each = N))
result <- adonis2(dist ~ group)
print(result$`Pr(>F)`[1])

# PERMANOVA, 2D
dist2 <- dist(pcoa$points[,1:2], method='euclidean')
result2 <- adonis2(dist2 ~ group)
print(result2$`Pr(>F)`[1])


# Save results
saveRDS(mat, sprintf("sim_rev_3-N%d-data.Rds", 2*N))
write.csv(group, sprintf('sim_rev-N%d-Y.csv', 2*N), row.names=FALSE)
saveRDS(dist, "sim_rev_1-dist.Rds")