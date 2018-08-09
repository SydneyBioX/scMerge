
findHVG <- function(exprsMat_HVG, batch,  intersection = 1, fdr = 0.01, minBiolDisp = 0.5){
  batch_list <- as.list(as.character(unique(batch)))
  HVG_list <- lapply(batch_list, function(x) {
    zeros <- apply(exprsMat_HVG[, batch == x], 1, function(x) sum(x == 0) / length(x))
    express_gene <- names(which(zeros <= 0.9))
    M3Drop::BrenneckeGetVariableGenes(exprsMat_HVG[express_gene, batch == x],
                                      suppress.plot = TRUE,
                                      fdr = 0.01,
                                      minBiolDisp = 0.5
    )
  })
  names(HVG_list) <- batch_list
  res <- unlist(HVG_list)
  tab <- table(res)
  HVG <- names(tab)[tab>=intersection]
  print(length(HVG))
  return(list(HVG = HVG, HVG_list = HVG_list))
}


identifyCluster <- function(exprsMat, batch, marker=NULL, HVG_list, kmeansK){

  batch_list <- as.list(as.character(unique(batch)))
  batch_oneType <- unlist(batch_list)[which(kmeansK == 1)]
  batch_num <- table(batch)[as.character(unique(batch))]

  if (!is.null(marker)) {
    pca <- lapply(
      batch_list,
      function(x) {
        if (!x %in% batch_oneType) {
          rsvd::rpca(t(exprsMat[marker, batch == x]), scale = T)
        } else {
          NULL
        }
      }
    )
  } else {
    pca <- lapply(
      batch_list,
      function(x) {
        if (!x %in% batch_oneType) {
          rsvd::rpca(t(exprsMat[HVG_list[[x]], batch == x]), scale = T)
        } else {
          NULL
        }
      }
    )
  }

  names(pca) <- unlist(batch_list)

  clustering_res <- list()
  clustering_res_pt_dist <- list()
  for (j in 1:length(pca)) {
    pca_current <- pca[[j]]
    if (!is.null(pca_current)) {
      kmeans_res <- kmeans(pca_current$x[, 1:10], centers = kmeansK[j], nstart = 100)
      clustering_res[[j]] <- kmeans_res$cluster
      clustering_res_pt_dist[[j]] <- lapply(1:kmeansK[j], function(y) {
        point_dist <- rowSums((pca[[j]]$x[which(kmeans_res$cluster == y), 1:10, drop = FALSE] - kmeans_res$centers[y, ])^2)
        point_rank <- rank(point_dist)
        point_rank <- point_rank / length(point_rank)
        point_rank
      })
      clustering_res_pt_dist[[j]] <- unlist(clustering_res_pt_dist[[j]])
      clustering_res_pt_dist[[j]] <- clustering_res_pt_dist[[j]][names(clustering_res[[j]])]
    } else {
      clustering_res[[j]] <- rep(1, batch_num[j])
      names(clustering_res[[j]]) <- colnames(exprsMat[, batch == batch_list[j]])
      # exprs_batch <- exprsMat[, batch == batch_list[j]]
      # centroid_batch <- rowMedians(exprs_batch)
      # point_dist <- colSums((exprs_batch - centroid_batch)^2)
      # point_rank <- rank(point_dist)
      # point_rank <- point_rank / length(point_rank)
      clustering_res_pt_dist[[j]] <- centroidDist(exprsMat[, batch == batch_list[j]])
    }
  }

  return(list(clustering_list = clustering_res, clustering_distProp = clustering_res_pt_dist))
}

centroidDist <- function(exprsMat){
  centroid_batch <- rowMedians(exprsMat)
  point_dist <- colSums((exprsMat - centroid_batch)^2)
  point_rank <- rank(point_dist)
  point_rank <- point_rank / length(point_rank)
  return(point_rank)
}

