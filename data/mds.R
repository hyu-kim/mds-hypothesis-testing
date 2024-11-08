data_df <- as.matrix(read.csv(sprintf('result/HyperparameterStudy/sim_%d-data.csv',1)))
y_df <- as.matrix(read.csv(sprintf('result/HyperparameterStudy/sim_%d-Y.csv',1)))
dist_mat <- as.matrix(dist(data_df[,1:4]))
z0 <- cmdscale(dist_mat, k = 2)

write.csv(z0, sprintf('result/Evaluation/sim_1-mds-Z.csv'), row.names=FALSE)
