# @created on: 10 July, 2018
# @last edited on: 10 July, 2018
# @author: Yingxin Lin

# require(ruv)
# require(igraph)
# require(pdist)
# require(proxy)


f_measure <- function(celltypes, batch) {
  f <- 2 * (celltypes * batch) / (celltypes + batch)
  return(f)
}

# Function to find mutual nearest clusters (for only two batches)
find_mnc <- function(sce_combine, res1, res2, exprs = "log2cpm", return_res = F, dist = "euclidean") {
  dist_res <- list()
  print("Calculating the distance")
  for (i in 1:max(res1)) {
    res_tmp <- c()
    for (j in 1:max(res2)) {
      mat1 <- assay(sce_combine[, names(which(res1 == i))], exprs)
      mat2 <- assay(sce_combine[, names(which(res2 == j))], exprs)
      if (dist == "cosine") {
        dist_mat <- dist(t(mat1), t(mat2), method = "cosine")
      } else if (dist == "cor") {
        dist_mat <- 1 - cor((mat1), (mat2))
      } else {
        dist_mat <- pdist(t(mat1), t(mat2))
      }
      dist_mat <- as.matrix(dist_mat)
      res_tmp <- c(res_tmp, mean(dist_mat))
    }
    dist_res[[i]] <- res_tmp
  }
  dist_res <- do.call(rbind, dist_res)
  # Find mutual nearest clusters
  print("Find mutual nearest clusters")
  neighbour_batch1 <- apply(dist_res, 1, which.min)
  neighbour_batch2 <- apply(dist_res, 2, which.min)
  mnc <- c()
  for (k in 1:length(neighbour_batch1)) {
    if (neighbour_batch2[neighbour_batch1[k]] == k) {
      mnc <- rbind(mnc, c(k, neighbour_batch1[k]))
    }
  }
  colnames(mnc) <- c("Batch1", "Batch2")
  if (return_res) {
    res <- c()
    res$neighbour_batch1 <- neighbour_batch1
    res$neighbour_batch2 <- neighbour_batch2
    res$mnc <- mnc
    return(res)
  } else {
    return(mnc)
  }
}


# Function to construct replicate matrix based on mnc results (for only two batches)

require(igraph)

mnc_replicate <- function(res1, res2, mnc) {
  if (!is.null(mnc)) {
    tmp <- matrix(0L,
      nrow = length(res1) + length(res2),
      ncol = max(res1) + max(res2) - nrow(mnc)
    )

    rownames(tmp) <- c(names(res1), names(res2))
    colnames(tmp) <- paste("replicate", 1:ncol(tmp), sep = "_")
    for (i in 1:nrow(mnc)) {
      tmp[c(names(which(res1 == mnc[i, 1])), names(which(res2 == mnc[i, 2]))), i] <- 1L
    }
    idx1_noRep <- which(!(1:max(res1)) %in% mnc[, 1])
    if (length(idx1_noRep) != 0) {
      for (j in 1:length(idx1_noRep)) {
        tmp[names(which(res1 == idx1_noRep[j])), nrow(mnc) + j] <- 1L
      }
    }

    idx2_noRep <- which(!(1:max(res2)) %in% mnc[, 2])
    if (length(idx2_noRep) != 0) {
      for (k in 1:length(idx2_noRep)) {
        tmp[names(which(res2 == idx2_noRep[k])), nrow(mnc) + length(idx1_noRep) + k] <- 1L
      }
    }
  } else {
    tmp <- matrix(0L,
      nrow = length(res1) + length(res2),
      ncol = max(res1) + max(res2)
    )
    rownames(tmp) <- c(names(res1), names(res2))
    colnames(tmp) <- paste("replicate", 1:ncol(tmp), sep = "_")

    for (j in 1:max(res1)) {
      tmp[names(which(res1 == j)), j] <- 1L
    }

    for (k in 1:max(res2)) {
      tmp[names(which(res2 == k)), max(res1) + k] <- 1L
    }
  }


  return(tmp)
}




# An extreamly complicated and unefficient function to find mnc for datasets with multiple batches
# mnc is combined by using network clustering
#
find_mnc_network <- function(sce_combine, kmeans_res_list, exprs = "log2cpm", dist = "euclidean") {
  # require(pdist)
  batch_num <- length(kmeans_res_list)
  names(kmeans_res_list) <- paste("Batch", 1:batch_num, sep = "")

  # Check whether there are batches that has only one cluster
  batch_oneType <- which(unlist(lapply(kmeans_res_list, function(x) length(levels(as.factor(x))) == 1)))

  if (length(batch_oneType) != 0) {
    if (length(batch_oneType) == batch_num) {
      combine_pair <- combn(batch_num, 2)
      batch_oneType <- NULL
      allones <- T
    } else {
      combine_pair <- combn(c(1:batch_num)[-batch_oneType], 2)
      for (i in batch_oneType) {
        for (j in c(1:batch_num)[-batch_oneType]) {
          combine_pair <- cbind(combine_pair, c(i, j))
        }
      }
      allones <- F
    }
  } else {
    combine_pair <- combn(batch_num, 2)
    allones <- F
  }
  # combine_pair <- combn(batch_num,2)




  mnc <- list()

  if (allones) {
    dist_res <- matrix(NA, nrow = batch_num, ncol = batch_num)
    for (k in 1:ncol(combine_pair)) {
      print(k)
      res1 <- kmeans_res_list[[combine_pair[1, k]]]
      res2 <- kmeans_res_list[[combine_pair[2, k]]]
      mat1 <- assay(sce_combine[, names(which(res1 == 1))], exprs)
      mat2 <- assay(sce_combine[, names(which(res2 == 1))], exprs)
      if (dist == "cosine") {
        # require(proxy)
        dist_mat <- dist(t(mat1), t(mat2), method = "cosine")
      } else if (dist == "cor") {
        dist_mat <- 1 - cor((mat1), (mat2))
      } else {
        dist_mat <- pdist(t(mat1), t(mat2))
      }
      dist_mat <- as.matrix(dist_mat)
      dist_res[combine_pair[1, k], combine_pair[2, k]] <- dist_res[combine_pair[2, k], combine_pair[1, k]] <- median(dist_mat)
    }
    neighbour_res <- apply(dist_res, 1, which.min)
    mnc_mat <- c()
    for (i in 1:length(neighbour_res)) {
      if (neighbour_res[neighbour_res[i]] == i) {
        mnc_mat <- rbind(mnc_mat, sort(c(i, neighbour_res[i])))
      }
    }
    mnc_mat <- unique(mnc_mat)
    mnc <- list()
    for (i in 1:nrow(mnc_mat)) {
      mnc[[i]] <- matrix(1, ncol = 2, nrow = 1)
      colnames(mnc[[i]]) <- c(
        paste("Batch", mnc_mat[i, 1], sep = ""),
        paste("Batch", mnc_mat[i, 2], sep = "")
      )
    }
  } else {
    for (k in 1:ncol(combine_pair)) {
      dist_res <- list()
      print(k)
      res1 <- kmeans_res_list[[combine_pair[1, k]]]
      res2 <- kmeans_res_list[[combine_pair[2, k]]]
      for (i in 1:max(res1)) {
        res_tmp <- c()
        for (j in 1:max(res2)) {
          mat1 <- assay(sce_combine[, names(which(res1 == i))], exprs)
          mat2 <- assay(sce_combine[, names(which(res2 == j))], exprs)
          if (dist == "cosine") {
            # require(proxy)
            dist_mat <- dist(t(mat1), t(mat2), method = "cosine")
          } else if (dist == "cor") {
            dist_mat <- 1 - cor((mat1), (mat2))
          } else {
            dist_mat <- pdist(t(mat1), t(mat2))
          }
          dist_mat <- as.matrix(dist_mat)
          res_tmp <- c(res_tmp, median(dist_mat))
        }
        dist_res[[i]] <- res_tmp
      }
      dist_res <- do.call(rbind, dist_res)
      neighbour_batch1 <- apply(dist_res, 1, which.min)
      neighbour_batch2 <- apply(dist_res, 2, which.min)
      mnc_tmp <- c()
      for (l in 1:length(neighbour_batch1)) {
        if (neighbour_batch2[neighbour_batch1[l]] == l) {
          mnc_tmp <- rbind(mnc_tmp, c(l, neighbour_batch1[l]))
        }
      }
      mnc[[k]] <- mnc_tmp
      colnames(mnc[[k]]) <- c(
        paste("Batch", combine_pair[1, k], sep = ""),
        paste("Batch", combine_pair[2, k], sep = "")
      )
    }
  }


  ### NEED FIX: NO MNC

  # Perform network analysis


  edge_list <- do.call(rbind, lapply(mnc, function(x)
    t(apply(x, 1, function(y) paste(colnames(x), y, sep = "_")))))


  g <- igraph::graph_from_edgelist(edge_list, directed = FALSE)
  plot(g)
  kc <- igraph::fastgreedy.community(g)
  kc_df <- data.frame(group = as.numeric(kc$membership), batch = as.numeric(gsub("Batch", "", gsub("_.*", "", kc$names))), cluster = as.numeric(gsub(".*_", "", kc$names)))

  if (allones) {
    kc_df_new <- kc_df
    batch_oneType <- c(1:batch_num)[-c(mnc_mat)]
    for (i in batch_oneType) {
      print(i)
      neighbour_order <- rank(dist_res[i, ], na.last = T)
      group_order1 <- kc_df[kc_df[, "batch"] == which(neighbour_order == 1), "group"]
      group_order2 <- kc_df[kc_df[, "batch"] == which(neighbour_order == 2), "group"]
      if (group_order1 == group_order2) {
        kc_df_new <- rbind(kc_df_new, c(group_order1, i, 1))
      }
    }
    kc_df <- kc_df_new
  }

  return(kc_df)
}

# An extreamly complicated and unefficient function to construct replicate matrix
# based on results from find_mnc_network




mnc_replicate_multiple_new <- function(kmeans_res_list, kc_df) {
  batch_num <- length(kmeans_res_list)
  if (!is.null(kc_df)) {
    idx_noRep <- list()

    # For each batches, check whether there is any clusters that do not have replicates
    # for(j in 1:length(kmeans_res_list)){
    #   idx_noRep[[j]]<-which(!(1:max(kmeans_res_list[[j]]))%in%mnc[,j])
    # }
    for (j in 1:length(kmeans_res_list)) {
      idx_noRep[[j]] <- which(!(1:max(kmeans_res_list[[j]])) %in% kc_df[kc_df$batch == j, "cluster"])
    }


    tmp <- matrix(0L,
      nrow = do.call(sum, lapply(kmeans_res_list, length)),
      ncol = max(kc_df$group) + do.call(sum, lapply(idx_noRep, length))
    )

    rownames(tmp) <- unlist(lapply(kmeans_res_list, names))
    colnames(tmp) <- paste("replicate", 1:ncol(tmp), sep = "_")

    # For each replicates
    replicate_size <- table(kc_df$group)
    for (i in 1:max(kc_df$group)) {
      tmp_names <- c()

      # For each batch in this replicate
      kc_df_sub <- kc_df[kc_df$group == i, ]
      for (l in 1:replicate_size[i]) {
        # tmp_names<-c(tmp_names,names(which(kmeans_res_list[[l]]==mnc[i,l])))

        tmp_names <- c(tmp_names, names(which(kmeans_res_list[[kc_df_sub[l, "batch"]]] == kc_df_sub[l, "cluster"])))
      }
      tmp[tmp_names, i] <- 1L
    }
    current_idx <- max(kc_df$group)
    for (j in 1:length(kmeans_res_list)) {
      if (length(idx_noRep[[j]]) != 0) {
        for (k in 1:length(idx_noRep[[j]])) {
          tmp[names(which(kmeans_res_list[[j]] == idx_noRep[[j]][k])), current_idx + k] <- 1L
        }
        current_idx <- current_idx + length(idx_noRep[[j]])
      }
    }
  }
  return(tmp)
}



mnc_replicate_multiple_new2 <- function(kmeans_res_list, kmeans_res_pt, replicate_prop, kc_df) {
  batch_num <- length(kmeans_res_list)
  if (!is.null(kc_df)) {
    idx_noRep <- list()

    # For each batches, check whether there is any clusters that do not have replicates
    # for(j in 1:length(kmeans_res_list)){
    #   idx_noRep[[j]]<-which(!(1:max(kmeans_res_list[[j]]))%in%mnc[,j])
    # }
    for (j in 1:length(kmeans_res_list)) {
      idx_noRep[[j]] <- which(!(1:max(kmeans_res_list[[j]])) %in% kc_df[kc_df$batch == j, "cluster"])
    }

    replicate_vector <- rep(NA, length(unlist(kmeans_res_list)))
    names(replicate_vector) <- names(unlist(kmeans_res_list))
    kmeans_res_pt <- unlist(kmeans_res_pt)
    replicate_size <- table(kc_df$group)
    for (i in 1:max(kc_df$group)) {
      tmp_names <- c()

      # For each batch in this replicate
      kc_df_sub <- kc_df[kc_df$group == i, ]
      for (l in 1:replicate_size[i]) {
        # tmp_names<-c(tmp_names,names(which(kmeans_res_list[[l]]==mnc[i,l])))

        tmp_names <- c(tmp_names, names(which(kmeans_res_list[[kc_df_sub[l, "batch"]]] == kc_df_sub[l, "cluster"])))
      }

      replicate_vector[tmp_names[kmeans_res_pt[tmp_names] <= replicate_prop]] <- paste("Replicate", i, sep = "_")
    }

    current_idx <- max(kc_df$group)
    for (j in 1:length(kmeans_res_list)) {
      if (length(idx_noRep[[j]]) != 0) {
        for (k in 1:length(idx_noRep[[j]])) {
          tmp_names <- names(which(kmeans_res_list[[j]] == idx_noRep[[j]][k]))
          replicate_vector[tmp_names[kmeans_res_pt[tmp_names] <= replicate_prop]] <- paste("Replicate", current_idx + k, sep = "_")
        }
        current_idx <- current_idx + length(idx_noRep[[j]])
      }
    }



    replicate_vector[is.na(replicate_vector)] <- 1:sum(is.na(replicate_vector))


    #
    #     tmp = matrix(0L, nrow = do.call(sum,lapply(kmeans_res_list,length)),
    #                  ncol = max(kc_df$group)+do.call(sum,lapply(idx_noRep,length)))
    #
    #     rownames(tmp) = unlist(lapply(kmeans_res_list,names))
    #     colnames(tmp) = paste("replicate",1:ncol(tmp),sep="_")
    #
    #     #For each replicates
    #     replicate_size<-table(kc_df$group)
    #     for(i in 1:max(kc_df$group)){
    #       tmp_names<-c()
    #
    #       #For each batch in this replicate
    #       kc_df_sub<-kc_df[kc_df$group==i,]
    #       for(l in 1:replicate_size[i]){
    #         # tmp_names<-c(tmp_names,names(which(kmeans_res_list[[l]]==mnc[i,l])))
    #
    #         tmp_names<-c(tmp_names,names(which(kmeans_res_list[[kc_df_sub[l,"batch"]]]==kc_df_sub[l,"cluster"])))
    #       }
    #       tmp[tmp_names,i]<-1L
    #     }
    #     current_idx<-max(kc_df$group)
    #     for(j in 1:length(kmeans_res_list)){
    #
    #       if(length(idx_noRep[[j]])!=0){
    #         for(k in 1:length(idx_noRep[[j]])){
    #           tmp[names(which(kmeans_res_list[[j]]==idx_noRep[[j]][k])),current_idx+k]<-1L
    #         }
    #         current_idx<-current_idx+length(idx_noRep[[j]])
    #       }
    #
    #     }
  }
  return(replicate_vector)
}


minMaxScale <- function(v) {
  v <- (v - min(v)) / (max(v) - min(v))
  return(v)
}

standardize <- function(dat, batch) {
  num_cell <- ncol(dat)
  batch <- as.factor(batch)
  grand.mean <- matrix(rowMeans(dat), nrow = 1)
  stand.mean <- t(grand.mean) %*% t(rep(1, num_cell))
  design <- model.matrix(~-1 + batch)
  B.hat <- solve(crossprod(design), tcrossprod(t(design), as.matrix(dat)))
  var.pooled <- ((dat - t(design %*% B.hat))^2) %*% rep(1 / num_cell, num_cell)
  s.data <- (dat - stand.mean) / (sqrt(var.pooled) %*% t(rep(1, num_cell)))
  return(res = list(s.data = s.data, stand.mean = stand.mean, stand.var = var.pooled))
}



# A function to perform loacation/scale adjustement to data as the input of RUVIII
# which also provides the option to select optimal RUVk according to the silhouette coefficient
normalise_RUVIII <- function(Y = Y, M = M, ctl = ctl,
                             fullalpha = NULL,
                             k = k, return.info = T,
                             cell_type = NULL,
                             batch = NULL,
                             return_all = F,
                             fast_svd = F) {
  Y <- t(Y)
  # geneMeanMat<-matrix(rep(rowMeans(Y),ncol(Y)),ncol=ncol(Y))
  # geneSdMat<-matrix(rep(apply(Y-geneMeanMat,1,sd),ncol(Y)),ncol=ncol(Y))
  # normY<-(Y-geneMeanMat)/geneSdMat

  if (is.null(batch)) {
    warning("No batch info!")
    return(NULL)
  }

  scale_res <- standardize(Y, batch)
  normY <- scale_res$s.data
  geneSdMat <- sqrt(scale_res$stand.var) %*% t(rep(1, ncol(Y)))
  geneMeanMat <- scale_res$stand.mean

  if (length(k) == 1) {
    if (fast_svd) {
      ruv3res <- fastRUVIII(
        Y = t(normY),
        ctl = ctl,
        k = k,
        M = M,
        fullalpha = fullalpha,
        return.info = return.info
      )
    } else {
      ruv3res <- RUVIII(
        Y = t(normY),
        ctl = ctl,
        k = k,
        M = M,
        fullalpha = fullalpha,
        return.info = return.info
      )
    }

    ruv3res$k <- k
  } else {
    ruv3res_list <- list()
    print(paste("k =", k[1]))
    if (fast_svd) {
      ruv3res_list[[1]] <- fastRUVIII(
        Y = t(normY),
        ctl = ctl,
        k = k[1],
        M = M,
        fullalpha = fullalpha,
        return.info = return.info
      )
    } else {
      ruv3res_list[[1]] <- RUVIII(
        Y = t(normY),
        ctl = ctl,
        k = k[1],
        M = M,
        fullalpha = fullalpha,
        return.info = return.info
      )
    }

    for (i in 2:length(k)) {
      print(paste("k =", k[i]))
      if (fast_svd) {
        ruv3res_list[[i]] <- fastRUVIII(
          Y = t(normY),
          ctl = ctl,
          k = k[i],
          M = M,
          fullalpha = ruv3res_list[[1]]$fullalpha,
          return.info = return.info
        )
      } else {
        ruv3res_list[[i]] <- RUVIII(
          Y = t(normY),
          ctl = ctl,
          k = k[i],
          M = M,
          fullalpha = ruv3res_list[[1]]$fullalpha,
          return.info = return.info
        )
      }
    }

    if (is.null(cell_type)) {
      cat("No cell type info, replicate matrix will be used as cell type info\n")
      cell_type <- apply(M, 1, function(x) which(x == 1))
    }

    sil_res <- do.call(cbind, lapply(ruv3res_list,
      FUN = function(x) {
        pca.data <- rpca(x$newY, k = 50, rand = 1)
        # pca.data<-prcomp(x$newY)
        c(
          batch_sil(pca.data, as.numeric(as.factor(cell_type))),
          batch_sil(pca.data, as.numeric(as.factor(batch)), nPCs = 10)
        )
      }
    ))
    # score<-minMaxScale(k*sil_res[1,])-minMaxScale(rev(k)*sil_res[2,])
    # f_score<-minMaxScale(sil_res[1,])-minMaxScale(sil_res[2,])
    f_score <- rep(NA, ncol(sil_res))
    for (i in 1:length(k)) {
      f_score[i] <- f_measure(minMaxScale(sil_res[1, ])[i], 1 - minMaxScale(sil_res[2, ])[i])
    }
    names(f_score) <- k


    ruv3res <- ruv3res_list[[which.max(f_score)]]
    ruv3res$k <- which.max(f_score)
    print(sil_res)
    print(paste("optimal k:", ruv3res$k))

    plot(k, f_score, pch = 16, col = "light grey")
    lines(k, f_score)
    points(ruv3res$k, f_score[ruv3res$k], col = "red", pch = 16)
  }
  if (return_all) {
    for (i in 1:length(ruv3res_list)) {
      ruv3res_list[[i]]$newY <- t((t(ruv3res_list[[i]]$newY) * geneSdMat + geneMeanMat))
    }
    if (is.null(batch)) {
      return(ruv3res_list)
    } else {
      ruv3res_list$optimise_k <- which.max(f_score)
      return(ruv3res_list)
    }
  } else {
    ruv3res$newY <- t((t(ruv3res$newY) * geneSdMat + geneMeanMat))
    return(ruv3res)
  }
}


wanted_variation_replicate <- function(sce, WV, WV_marker, replicate_vector, exprs) {
  names(WV) <- colnames(sce)
  if (!is.null(WV_marker)) {
    marker_expr <- lapply(WV_marker, function(x) aggregate(assay(sce, exprs)[x, ],
        by = list(replicate_vector),
        FUN = mean
      ))
    names(marker_expr) <- WV_marker

    marker_expr <- do.call(cbind, lapply(marker_expr, function(x) x[, 2]))
    rownames(marker_expr) <- names(table(replicate_vector))
    km <- kmeans(marker_expr, centers = 2)
    tab_max <- table(apply(km$centers, 2, which.max))
    marker_cluster <- names(tab_max)[which.max(tab_max)]
    marker_replicate <- names(which(km$cluster == marker_cluster))
    replicate_vector[replicate_vector %in% marker_replicate] <- paste(replicate_vector[replicate_vector %in% marker_replicate],
      WV[replicate_vector %in% marker_replicate],
      sep = "_"
    )
    repMat_new <- ruv::replicate.matrix(as.factor(replicate_vector))
  } else {
    repMat_new <- ruv::replicate.matrix(as.factor(paste(replicate_vector, WV, sep = "_")))
  }

  return(repMat_new)
}



# if (!is.null(cell_type) & is.null(cell_type_inc) & !cell_type_match) {
#
#   cat("Performing supervised scMerge\n")
#   names(cell_type) <- colnames(exprs_mat)
#
#   repVector <- supervisedReplicate(exprs_mat, cell_type, replicate_prop)
#
#   repMat <- ruv::replicate.matrix(repVector)
#
# } else if (!is.null(cell_type) & is.null(cell_type_inc) & cell_type_match) {
#
#   cat("Finding MNC from the known cell types of different batches...\n")
#   names(cell_type) <- colnames(exprs_mat)
#   sce_batch <- as.factor(sce_batch)
#   batch_list <- as.list(as.character(unique(sce_batch)))
#
#
#   if (is.null(marker)) {
#
#     cat("Finding HVG...\n")
#     exprsMat_HVG <- assay(sce_combine, hvg_exprs)
#     HVG_res <- findHVG(exprsMat_HVG, sce_batch)
#     Brennecke_HVG <- HVG_res$HVG
#     Brennecke_HVG_list <- HVG_res$HVG_list
#
#   } else {
#     Brennecke_HVG <- marker
#   }
#   # sce_combine$batch <- as.factor(sce_combine$batch)
#   # batch_list <- as.list(as.character(unique(sce_combine$batch)))
#   clustering_distProp <- lapply(unique(cell_type),
#                                 function(x) centroidDist(exprs_mat[,cell_type==x]))
#   clustering_distProp <- unlist(clustering_distProp)[names(cell_type)]
#
#   clustering_distProp_list_batch <- lapply(batch_list, function(x) {
#     tmp <- clustering_distProp[sce_batch == x]
#     names(tmp) <- colnames(exprs_mat)[sce_batch == x]
#     tmp
#   })
#
#
#   cellType_list_batch <- lapply(batch_list, function(x) {
#     tmp <- as.numeric(droplevels(as.factor(cell_type[sce_batch == x])))
#     names(tmp) <- colnames(exprs_mat)[sce_batch == x]
#     tmp
#   })
#
#   mnc_res <- findMNC(exprs_mat[Brennecke_HVG, ],
#                      clustering_list = cellType_list_batch,
#                      dist = dist)
#   print(mnc_res)
#
#   repVector <- mncRepcliate(clustering_list = cellType_list_batch,
#                             clustering_distProp = clustering_distProp_list_batch,
#                             replicate_prop = replicate_prop,
#                             mnc_df = mnc_res)
#   repMat <- ruv::replicate.matrix(repVector)
#
#
#   # if(!is.null(cell_type)&!is.null(cell_type_inc)){
#   #   cat("Performing semi-supervised scMerge with subsets of known cell type\n")
#   #   repVector[cell_type_inc] <- cell_type[cell_type_inc]
#   # }
#   # repMat <- replicate.matrix(repVector)
# } else {
#
#   sce_batch <- as.factor(sce_batch)
#   batch_list <- as.list(as.character(unique(sce_batch)))
#
#
#   if (is.null(kmeansK)) {
#     stop("KmeansK is NULL", call. = FALSE)
#   }
#
#   if (length(batch_list) != length(kmeansK)) {
#     stop("length of KmeansK needs to be the same as the number of batch", call. = FALSE)
#   }
#
#
#   # Find HVG
#   if (is.null(marker)) {
#     cat("Finding HVG...\n")
#     exprsMat_HVG <- assay(sce_combine, hvg_exprs)
#     HVG_res <- findHVG(exprsMat_HVG, sce_batch)
#     Brennecke_HVG <- HVG_res$HVG
#     Brennecke_HVG_list <- HVG_res$HVG_list
#   } else {
#     Brennecke_HVG <- marker
#   }
#
#
#   # Clustering within each batch
#
#   cat("Clustering within each batch...\n")
#
#   cat("Performing pca...\n")
#
#   cluster_res <- identifyCluster(exprsMat = exprs_mat,
#                                  batch = sce_batch,
#                                  marker = marker,
#                                  HVG_list = Brennecke_HVG_list,
#                                  kmeansK = kmeansK)
#
#   # Find Mutual Nearest Cluster
#
#   cat("Creating Mutual Nearest Cluster...\n")
#
#   mnc_res <- findMNC(exprs_mat[Brennecke_HVG, ],
#                      clustering_list = cluster_res$clustering_list,
#                      dist = dist)
#
#   print(mnc_res)
#
#   # Create replicate matrix
#   repVector <- mncRepcliate(clustering_list = cluster_res$clustering_list,
#                             clustering_distProp = cluster_res$clustering_distProp,
#                             replicate_prop = replicate_prop,
#                             mnc_df = mnc_res)
#
#   if (!is.null(cell_type) & !is.null(cell_type_inc)) {
#     cat("Performing semi-supervised scMerge with subsets of known cell type\n")
#     repVector[cell_type_inc] <- cell_type[cell_type_inc]
#   }
#   repMat <- ruv::replicate.matrix(repVector)
#
#   if (!is.null(WV)) {
#     cat("Performing semi-supervised scMerge with wanted variation\n")
#     repMat <- wvReplicate(exprs_mat, WV, WV_marker, repVector)
#   }
# }
