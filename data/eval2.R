rp_eval2 <- function(nperm, dm, y_orig, emb, y_emb = NULL){
  n <- length(y_orig)
  if(is.null(y_emb)){
    y_emb <- y_orig
  }
  F0 <- pseudo_F(trt = y_orig, d = dm)$ratio
  Femb <- pseudo_F(mat = emb, trt = y_emb)$ratio
  
  F0_perm <- Femb_perm <- rep(0, nperm)
  for(i in 1:nperm){
    set.seed(i)
    ind_perm <- sample(1:n, size = n, replace = F)
    y0_perm <- y_orig[ind_perm]
    ye_perm <- y_emb[ind_perm]
    
    F0_perm[i] <- pseudo_F(trt = y0_perm, d = dm)$ratio
    Femb_perm[i] <- pseudo_F(mat = emb, trt = ye_perm)$ratio
  }
  
  q = which(Femb == sort(c(Femb_perm, Femb))) /
    which(F0 == sort(c(F0_perm, F0)))
  return(list(F0 = F0, F0_perm = F0_perm,
              Femb = Femb, Femb_perm = Femb_perm,
              rho = cor(F0 - F0_perm, Femb - Femb_perm),
              q = q))
}

### Eval
## Cirrhosis
cirr_res$mds <- data.frame(cmdscale(phyl_unifrac_cirrhosis, k = 2))

eval_res_cirr <- list(
  mds = rp_eval2(500, phyl_unifrac_cirrhosis, cirr_y, cirr_res$mds),
  fmds3 = rp_eval2(500, phyl_unifrac_cirrhosis, cirr_y, cirr_res$fmds$lambda0.3$z),
  fmds5 = rp_eval2(500, phyl_unifrac_cirrhosis, cirr_y, cirr_res$fmds$lambda0.5$z), 
  fmds7 = rp_eval2(500, phyl_unifrac_cirrhosis, cirr_y, cirr_res$fmds$lambda0.7$z), 
  
  umap_s5 = rp_eval2(500, phyl_unifrac_cirrhosis, cirr_y, cirr_res$umap[[1]]), 
  umap_s10 = rp_eval2(500, phyl_unifrac_cirrhosis, cirr_y, cirr_res$umap[[2]]), 
  umap_s20 = rp_eval2(500, phyl_unifrac_cirrhosis, cirr_y, cirr_res$umap[[3]]), 
  umap_s30 = rp_eval2(500, phyl_unifrac_cirrhosis, cirr_y, cirr_res$umap[[4]]), 
  
  umap_u5 = rp_eval2(500, phyl_unifrac_cirrhosis, cirr_y, cirr_res$umap[[5]]), 
  umap_u10 = rp_eval2(500, phyl_unifrac_cirrhosis, cirr_y, cirr_res$umap[[6]]), 
  umap_u20 = rp_eval2(500, phyl_unifrac_cirrhosis, cirr_y, cirr_res$umap[[7]]), 
  umap_u30 = rp_eval2(500, phyl_unifrac_cirrhosis, cirr_y, cirr_res$umap[[8]]), 
  
  tsne5 =  rp_eval2(500, phyl_unifrac_cirrhosis, cirr_y, cirr_res$tsne[[5]]), 
  tsne7 =  rp_eval2(500, phyl_unifrac_cirrhosis, cirr_y, cirr_res$tsne[[7]]), 
  tsne10 =  rp_eval2(500, phyl_unifrac_cirrhosis, cirr_y, cirr_res$tsne[[10]]), 
  
  iso5 = rp_eval2(500, phyl_unifrac_cirrhosis, cirr_y, cirr_res$isomap[[5]]$points), 
  iso7 = rp_eval2(500, phyl_unifrac_cirrhosis, cirr_y, cirr_res$isomap[[7]]$points), 
  iso10 = rp_eval2(500, phyl_unifrac_cirrhosis, cirr_y, cirr_res$isomap[[10]]$points), 
)



# library(superMDS)
# sim_res$smds <- list(
#   TrainSuperMDS(d = sim_data$distmat, y = sim_data$data$Y, alpha = 0.1),
#   TrainSuperMDS(d = sim_data$distmat, y = sim_data$data$Y, alpha = 0.3),
#   TrainSuperMDS(d = sim_data$distmat, y = sim_data$data$Y, alpha = 0.5)
# )
# eval_res2$sim$smds <- NULL





par(mfrow = c(2, 3), mar = c(4,4,2,1))
for(i in c(1:23)){
  qqplot(eval_res2$sim[[i]]$F0_perm, eval_res2$sim[[i]]$Femb_perm, 
         main = names(eval_res2$sim)[i])
  abline(a=0, b=1, col = "red")
}

par(mfrow = c(2, 3), mar = c(4,4,2,1))
for(i in c(1:23)){
  plot(eval_res2$sim[[i]]$F0_perm, eval_res2$sim[[i]]$Femb_perm, 
         main = names(eval_res2$sim)[i])
}

####### Site 1
eval_res2$site1 <- list(
  mds =   rp_eval2(500, distmat1, y1, zmds1),
  fmds1 = rp_eval2(500, distmat1, y1, res1$lambda0.1$z),
  fmds3 = rp_eval2(500, distmat1, y1, res1$lambda0.3$z),
  fmds5 = rp_eval2(500, distmat1, y1, res1$lambda0.5$z), 

  umap_s5 =  rp_eval2(500, distmat1, y1, umap1[[1]]), 
  umap_s10 = rp_eval2(500, distmat1, y1, umap1[[2]]), 
  umap_s20 = rp_eval2(500, distmat1, y1, umap1[[3]]), 
  umap_s30 = rp_eval2(500, distmat1, y1, umap1[[4]]),
  
  umap_u5 =  rp_eval2(500, distmat1, y1, umap1[[5]]), 
  umap_u10 = rp_eval2(500, distmat1, y1, umap1[[6]]), 
  umap_u20 = rp_eval2(500, distmat1, y1, umap1[[7]]), 
  umap_u30 = rp_eval2(500, distmat1, y1, umap1[[8]]),
  
  tsne5 =  rp_eval2(500, distmat1, y1, tsne1[[5]]),
  tsne7 =  rp_eval2(500, distmat1, y1, tsne1[[7]]),
  tsne10 = rp_eval2(500, distmat1, y1, tsne1[[10]]), 
  
  iso5 =  rp_eval2(500, distmat1, y1, isomap1[[1]]$points), 
  iso7 =  rp_eval2(500, distmat1, y1, isomap1[[2]]$points), 
  iso10 = rp_eval2(500, distmat1, y1, isomap1[[3]]$points), 
  
  smds1 = rp_eval2(500, distmat1, y1, smds_site1$a.1$z),
  smds3 = rp_eval2(500, distmat1, y1, smds_site1$a.3$z),
  smds5 = rp_eval2(500, distmat1, y1, smds_site1$a.5$z)
)


par(mfrow = c(2, 3), mar = c(4,4,2,1))
for(i in c(1:23)){
  qqplot(eval_res2$site1[[i]]$F0_perm, eval_res2$sim[[i]]$Femb_perm,
         main = names(eval_res2$sim)[i])
  abline(a=0, b=1, col = "red")
}

par(mfrow = c(2, 3), mar = c(4,4,2,1))
for(i in c(1:23)){
  plot(eval_res2$site1[[i]]$F0_perm, eval_res2$site1[[i]]$Femb_perm, 
       main = names(eval_res2$site1)[i])
}


####### Site 2
eval_res2$site2 <- list(
  mds =   rp_eval2(500, distmat2, y2, zmds2),
  fmds1 = rp_eval2(500, distmat2, y2, res2$lambda0.1$z),
  fmds3 = rp_eval2(500, distmat2, y2, res2$lambda0.3$z),
  fmds5 = rp_eval2(500, distmat2, y2, res2$lambda0.5$z), 
  
  umap_s5 =  rp_eval2(500, distmat2, y2, umap2[[1]]), 
  umap_s10 = rp_eval2(500, distmat2, y2, umap2[[2]]), 
  umap_s20 = rp_eval2(500, distmat2, y2, umap2[[3]]), 
  umap_s30 = rp_eval2(500, distmat2, y2, umap2[[4]]),
  
  umap_u5 =  rp_eval2(500, distmat2, y2, umap2[[5]]), 
  umap_u10 = rp_eval2(500, distmat2, y2, umap2[[6]]), 
  umap_u20 = rp_eval2(500, distmat2, y2, umap2[[7]]), 
  umap_u30 = rp_eval2(500, distmat2, y2, umap2[[8]]),
  
  tsne5 =  rp_eval2(500, distmat2, y2, tsne2[[5]]),
  tsne7 =  rp_eval2(500, distmat2, y2, tsne2[[7]]),
  tsne10 = rp_eval2(500, distmat2, y2, tsne2[[10]]), 
  
  iso5 =  rp_eval2(500, distmat2, y2, isomap2[[1]]$points), 
  iso7 =  rp_eval2(500, distmat2, y2, isomap2[[2]]$points), 
  iso10 = rp_eval2(500, distmat2, y2, isomap2[[3]]$points), 
  
  smds1 = rp_eval2(500, distmat2, y2, smds_site2$a.1$z),
  smds3 = rp_eval2(500, distmat2, y2, smds_site2$a.3$z),
  smds5 = rp_eval2(500, distmat2, y2, smds_site2$a.5$z)
)

par(mfrow = c(2, 3), mar = c(4,4,2,1))
# for(i in c(1:23)){
#   qqplot(eval_res2$site1[[i]]$F0_perm, eval_res2$sim[[i]]$Femb_perm, 
#          main = names(eval_res2$sim)[i])
#   abline(a=0, b=1, col = "red")
# }

par(mfrow = c(2, 3), mar = c(4,4,2,1))
for(i in c(1:21)){
  plot(eval_res2$site2[[i]]$F0_perm, eval_res2$site2[[i]]$Femb_perm, 
       main = names(eval_res2$site2)[i])
}






sapply(eval_res2$sim,   function(x){round(x$rho, 4)})
sapply(eval_res2$site1, function(x){round(x$rho, 4)})
sapply(eval_res2$site2, function(x){round(x$rho, 4)})


#########etc
nn1 <- cbind(
  read.csv("data//result//nn//labels_out_site1.csv", header = F),
  read.csv("data//result//nn//features_out_site1.csv", header = F)
)
names(nn1)[1] <- c("Y")
nn1

nn2 <- cbind(
  read.csv("data//result//nn//labels_out_site2.csv", header = F),
  read.csv("data//result//nn//features_out_site2.csv", header = F)
)
names(nn2)[1] <- c("Y")
nn2

eval_res2$site1

eval_res2$site1$nn <- rp_eval2(500, distmat1, y1, nn1[, 2:33])
eval_res2$site2$nn <- rp_eval2(500, distmat2, y2, nn2[, 2:33])

### other measures\

site1_embed <- list(
  mds =   zmds1,
  fmds1 = res1$lambda0.1$z,
  fmds3 = res1$lambda0.3$z,
  fmds5 = res1$lambda0.5$z, 
  
  umap_s5 =  umap1[[1]], 
  umap_s10 = umap1[[2]], 
  umap_s20 = umap1[[3]], 
  umap_s30 = umap1[[4]],
  
  umap_u5 =  umap1[[5]], 
  umap_u10 = umap1[[6]], 
  umap_u20 = umap1[[7]], 
  umap_u30 = umap1[[8]],
  
  tsne5 =  tsne1[[5]],
  tsne7 =  tsne1[[7]],
  tsne10 = tsne1[[10]], 
  
  iso5 =  isomap1[[1]]$points, 
  iso7 =  isomap1[[2]]$points, 
  iso10 = isomap1[[3]]$points, 
  
  smds1 = smds_site1$a.1$z,
  smds3 = smds_site1$a.3$z,
  smds5 = smds_site1$a.5$z,
  
  nn = nn1[, 2:33]
)

site2_embed <- list(
  mds =   zmds2,
  fmds1 = res2$lambda0.1$z,
  fmds3 = res2$lambda0.3$z,
  fmds5 = res2$lambda0.5$z, 
  
  umap_s5 =  umap2[[1]], 
  umap_s10 = umap2[[2]], 
  umap_s20 = umap2[[3]], 
  umap_s30 = umap2[[4]],
  
  umap_u5 =  umap2[[5]], 
  umap_u10 = umap2[[6]], 
  umap_u20 = umap2[[7]], 
  umap_u30 = umap2[[8]],
  
  tsne5 =  tsne2[[5]],
  tsne7 =  tsne2[[7]],
  tsne10 = tsne2[[10]], 
  
  iso5 =  isomap2[[1]]$points, 
  iso7 =  isomap2[[2]]$points, 
  iso10 = isomap2[[3]]$points, 
  
  smds1 = smds_site2$a.1$z,
  smds3 = smds_site2$a.3$z,
  smds5 = smds_site2$a.5$z,
  
  nn = nn2[, 2:33]
)

eval_res2$site1 <- lapply(site1_embed, function(x){
  rp_eval2(500, distmat1, y1, x, y1)
})
eval_res2$site1$nn <- rp_eval2(500, distmat1, y1, nn1[,2:33], nn1$Y)

eval_res2$site2 <- lapply(site2_embed, function(x){
  rp_eval2(500, distmat2, y2, x, y1)
})
eval_res2$site2$nn <- rp_eval2(500, distmat2, y2, nn2[,2:33], nn2$Y)


t(rbind(
  data.frame(lapply(eval_res2$site1, function(x){x$q})),
  data.frame(lapply(eval_res2$site1, function(x){x$rho})),
  data.frame(lapply(site1_embed, function(x){get_stress(dist1, dist(x))})),
  data.frame(lapply(site1_embed, function(x){cor(dist1, dist(x))}))
))


t(rbind(
  data.frame(lapply(eval_res2$site2, function(x){x$q})),
  data.frame(lapply(eval_res2$site2, function(x){x$rho})),
  data.frame(lapply(site2_embed, function(x){get_stress(dist2, dist(x))})),
  data.frame(lapply(site2_embed, function(x){cor(dist2, dist(x))}))
))


t(rbind(
  data.frame(lapply(eval_res2$site2, function(x){x$rho})),
  data.frame(lapply(site2_embed, function(x){get_stress(dist2, dist(x))})),
  data.frame(lapply(site2_embed, function(x){cor(dist2, dist(x))}))
))


sim_embed <- list(
  mds =   sim_res$mds,
  fmds3 = sim_res$proposed$lambda0.3$z,
  fmds5 = sim_res$proposed$lambda0.5$z, 
  fmds7 = sim_res$proposed$lambda0.7$z, 
  
  umap_s5 =  sim_res$umap[[1]], 
  umap_s10 = sim_res$umap[[2]], 
  umap_s20 = sim_res$umap[[3]], 
  umap_s30 = sim_res$umap[[4]],
  umap_u5 =  sim_res$umap[[5]], 
  umap_u10 = sim_res$umap[[6]], 
  umap_u20 = sim_res$umap[[7]], 
  umap_u30 = sim_res$umap[[8]],
  
  tsne5 =  sim_res$tsne[[1]],
  tsne10 = sim_res$tsne[[2]], 
  tsne20 = sim_res$tsne[[3]], 
  tsne30 = sim_res$tsne[[4]],
  iso5 =   sim_res$isomap[[1]]$points, 
  iso10 =  sim_res$isomap[[2]]$points, 
  iso20 =  sim_res$isomap[[3]]$points, 
  iso30 =  sim_res$isomap[[4]]$points,

  smds1 = sim_res$smds[[1]]$z,
  smds3 = sim_res$smds[[2]]$z,
  smds5 = sim_res$smds[[3]]$z
)

t(rbind(
  data.frame(lapply(eval_res2$sim, function(x){x$q})),
  data.frame(lapply(eval_res2$sim, function(x){x$rho})),
  data.frame(lapply(sim_embed, function(x){get_stress(sim_data$dist, dist(x))})),
  data.frame(lapply(sim_embed, function(x){cor(sim_data$dist, dist(x))}))
))

par(mfrow = c(1,1))
plot(sim_data$dist, dist(sim_embed$tsne30))
