#' Create pseudo-replicate for scMerge
#'
#' Create pseudo-replicate from single-cell RNA-seq data from different batches, experiments, and protocols.
#'
#' @param sce A \code{SingleCellExperiment} object contains the batch-combined matrix with batch info in colData
#' @param batch A vector indicates the batch information for each cell in the batch-combined matrix.
#' @param kmeansK A vector indicates the kmeans's K for each batch, length of KmeansK needs to be the same as the number of batch.
#' @param exprs A string inciates the assay that are used for batch correction, default is logcounts
#' @param hvg_exprs A string inciates the assay that are used for highly variable genes identification, default is counts
#' @param marker A vector of markers, which will be used in calcualtion of mutual nearest cluster. If no markers input, highly variable genes will be used instead
#' @param marker_list A list of markers for each batch, which will be used in calcualtion of mutual nearest cluster.
#' @param replicate_prop A number indicates the ratio of cells that are included in pseudo-replicates, ranges from 0 to 1
#' @param cell_type A vector indicates the cell type information for each cell in the batch-combined matrix. If it is \code{NULL}, pseudo-replicate procedure will be run to identify cell type.
#' @param cell_type_match Whether find mutual nearest cluster using cell type information
#' @param cell_type_inc A vector indicates the indices of the cells that will be used to supervise the pseudo-replicate procedure
#' @param dist The distance metrics that are used in the calcualtion of the mutual nearest cluster, default is Pearson correlation.
#' @param WV A vector indicates the wanted variation factor other than cell type info, such as cell stages.
#' @param WV_marker A vector indicates the markers of the wanted variation.
#' @param parallelParam The \code{BiocParallelParam} class from the \code{BiocParallel} package is used. Default is bpparam(),
#' @param return_all If \code{FALSE}, only return the replicate matrix.
#' @import BiocParallel
#'
#' @return If \code{return_all} is \code{FALSE}, return a replicate matrix.
#' If \code{return_sce} is \code{TRUE}, return the followings
#' \item{repMat }{replicate matrix}
#' \item{mnc }{mutual nearest cluster}
#' \item{replicate vector }{replicate vector}
#' \item{HVG }{highly variable genes used in scReplicate}
#' @author Yingxin Lin, Kevin Wang
#' @export
#' @examples
#' \dontrun{
#' suppressPackageStartupMessages({
#' library(SingleCellExperiment)
#' library(scMerge)
#' library(scMerge.data)
#' })
#' ## Loading example data
#' data("sce_mESC", package = "scMerge.data")
#'
#' scRep_result <- scReplicate(
#'   sce = sce_mESC,
#'   batch = sce_mESC$batch,
#'   kmeansK = c(1,3,3,1,1),
#'   fast_svd = FALSE)
#'   }


scReplicate <- function(sce, batch = NULL, kmeansK = NULL, exprs = "logcounts", hvg_exprs = "counts",
                        marker = NULL, marker_list = NULL, replicate_prop = 1, cell_type = NULL, cell_type_match = FALSE,
                        cell_type_inc = NULL, dist = "cor", WV = NULL, WV_marker = NULL,
                        parallelParam = bpparam(), return_all = FALSE, fast_svd) {

  exprs_mat <- SummarizedExperiment::assay(sce, exprs)

  ## Based on the input cell_type info, we decide which version of scMerge to use.
  if (!is.null(cell_type) & is.null(cell_type_inc) & !cell_type_match) {
    cat("Cell type information was supplied, and mutual nearest neighbour option was not selected \n")
    cat("Performing supervised scMerge \n")
    names(cell_type) <- colnames(exprs_mat)

    repVector <- supervisedReplicate(exprs_mat, cell_type, replicate_prop)

    mnc_res <- NULL


  } else if (!is.null(cell_type) & is.null(cell_type_inc) & cell_type_match) {
    cat("Cell type information was supplied, and mutual nearest neighbour option was selected \n")
    cat("Finding MNC from the known cell types of different batches...\n")
    names(cell_type) <- colnames(exprs_mat)
    batch <- as.factor(batch)
    batch_list <- as.list(as.character(unique(batch)))

    ## Finding data-adaptive markers, store as object as HVG
    if (is.null(marker) & is.null(marker_list)) {
      cat("No maker nor marker_list information was supplied \n")
      cat("Finding HVG...\n")
      exprs_mat_HVG <- SummarizedExperiment::assay(sce, hvg_exprs)
      HVG_res <- findHVG(exprs_mat_HVG, batch, parallelParam = parallelParam)
      HVG <- HVG_res$HVG
      HVG_list <- HVG_res$HVG_list

    } else if (!is.null(marker_list) & is.null(marker)) {
      ## Taking the union of marker_list as HVG
      cat("Only marker_list information was supplied \n")
      HVG_list <- marker_list
      names(HVG_list) <- batch_list
      HVG <- Reduce(union, marker_list)
    } else {
      ## If there is only a list of markers, then use that as HVG
      HVG <- marker
    }
    # sce$batch <- as.factor(sce$batch) batch_list <-
    # as.list(as.character(unique(sce$batch)))

    ## For every unique cell_type, calculate the centroid.
    clustering_distProp <- lapply(unique(cell_type), function(x) centroidDist(exprs_mat[,
                                                                                        cell_type == x]))
    clustering_distProp <- unlist(clustering_distProp)[names(cell_type)]

    ## For each batch in batch_list, we find the corresponding cell_type_centroid
    clustering_distProp_list_batch <- lapply(batch_list, function(x) {
      tmp <- clustering_distProp[batch == x]
      names(tmp) <- colnames(exprs_mat)[batch == x]
      tmp
    })

    ## For each batch in batch_list, we find the corresponding cell_type
    cellType_list_batch <- lapply(batch_list, function(x) {
      tmp <- as.numeric(droplevels(as.factor(cell_type[batch == x])))
      names(tmp) <- colnames(exprs_mat)[batch == x]
      tmp
    })

    ## Here, using the cell type information, we go ahead and find MNC
    mnc_res <- findMNC(exprs_mat[HVG, ], clustering_list = cellType_list_batch,
                                 dist = dist, parallelParam = parallelParam)

    print(mnc_res)

    repVector <- mncReplicate(clustering_list = cellType_list_batch, clustering_distProp = clustering_distProp_list_batch,
                              replicate_prop = replicate_prop, mnc_df = mnc_res)



    # if(!is.null(cell_type)&!is.null(cell_type_inc)){ cat('Performing
    # semi-supervised scMerge with subsets of known cell type\n')
    # repVector[cell_type_inc] <- cell_type[cell_type_inc] } repMat <-
    # replicate.matrix(repVector)
  } else {

    batch <- as.factor(batch)
    batch_list <- as.list(as.character(unique(batch)))

    ## KEVIN: we should move this input check much earlier, into scMerge function

    if (is.null(kmeansK)) {
      stop("KmeansK is NULL", call. = FALSE)
    }

    if (length(batch_list) != length(kmeansK)) {
      stop("length of KmeansK needs to be the same as the number of batch",
           call. = FALSE)
    }


    ## If there are no marker information what so ever, we find the HVG adaptively
    ## from the data
    if (is.null(marker) & is.null(marker_list)) {
      cat("No maker nor marker_list information was supplied \n")
      cat("Finding HVG...\n")
      exprs_mat_HVG <- SummarizedExperiment::assay(sce, hvg_exprs)
      HVG_res <- findHVG(exprs_mat_HVG, batch, parallelParam = parallelParam)
      HVG <- HVG_res$HVG
      HVG_list <- HVG_res$HVG_list

    } else if (!is.null(marker_list) & is.null(marker)) {
      cat("Only marker_list information was supplied \n")
      HVG_list <- marker_list
      names(HVG_list) <- batch_list
      HVG <- Reduce(union, marker_list)
    } else {
      HVG <- marker
    }


    # Clustering within each batch

    cat("Clustering within each batch...\n")

    cat("Performing pca...\n")

    cluster_res <- identifyCluster(exprs_mat = exprs_mat, batch = batch, marker = marker,
                                   HVG_list = HVG_list, kmeansK = kmeansK, parallelParam = parallelParam, fast_svd = fast_svd)

    # Find Mutual Nearest Cluster

    cat("Creating Mutual Nearest Cluster...\n")

    # if (parallel) {
    #   mnc_res <- findMNC_parallel(exprs_mat[HVG, ], clustering_list = cluster_res$clustering_list,
    #                               dist = dist, parallelParam = parallelParam)
    # } else {
    #   mnc_res <- findMNC(exprs_mat[HVG, ], clustering_list = cluster_res$clustering_list,
    #                      dist = dist)
    # }

    mnc_res <- findMNC(exprs_mat[HVG, ], clustering_list = cluster_res$clustering_list,
                                dist = dist, parallelParam = parallelParam)


    print(mnc_res)

    # Create replicate matrix
    repVector <- mncReplicate(clustering_list = cluster_res$clustering_list,
                              clustering_distProp = cluster_res$clustering_distProp, replicate_prop = replicate_prop,
                              mnc_df = mnc_res)

    if (!is.null(cell_type) & !is.null(cell_type_inc)) {
      cat("Performing semi-supervised scMerge with subsets of known cell type\n")
      repVector[cell_type_inc] <- cell_type[cell_type_inc]
    }


    if (!is.null(WV)) {
      cat("Performing semi-supervised scMerge with wanted variation\n")
      repVector <- wvReplicate(exprs_mat, WV, WV_marker, repVector)
    }


  }

  repMat <- ruv::replicate.matrix(repVector)


  if (return_all) {
    return(list(repMat = repMat, mnc = mnc_res, repVector = repVector, HVG = HVG))
  } else {
    return(repMat)
  }
}

########################### Function to find HVG ######################################################################


findHVG <- function(exprs_mat_HVG, batch, intersection = 1, fdr = 0.01, minBiolDisp = 0.5,
                    parallelParam) {
  batch_list <- as.list(as.character(unique(batch)))
  # if (parallel) {
    HVG_list <- BiocParallel::bplapply(batch_list, function(x) {
      zeros <- apply(exprs_mat_HVG[, batch == x], 1, function(x) sum(x == 0)/length(x))
      express_gene <- names(which(zeros <= 0.9))
      M3Drop::BrenneckeGetVariableGenes(exprs_mat_HVG[express_gene, batch == x], suppress.plot = TRUE,
                                        fdr = 0.01, minBiolDisp = 0.5)
    }, BPPARAM = parallelParam)
  # } else {
  #   HVG_list <- lapply(batch_list, function(x) {
  #     zeros <- apply(exprs_mat_HVG[, batch == x], 1, function(x) sum(x == 0)/length(x))
  #     express_gene <- names(which(zeros <= 0.9))
  #     M3Drop::BrenneckeGetVariableGenes(exprs_mat_HVG[express_gene, batch == x], suppress.plot = TRUE,
  #                                       fdr = 0.01, minBiolDisp = 0.5)
  #   })
  # }

  names(HVG_list) <- batch_list
  res <- unlist(HVG_list)
  tab <- table(res)
  HVG <- names(tab)[tab >= intersection]
  print(length(HVG))
  return(list(HVG = HVG, HVG_list = HVG_list))
}

###################################################################################################### Function to identify clusters from each batch
identifyCluster <- function(exprs_mat, batch, marker = NULL, HVG_list, kmeansK, parallelParam,
                            fast_svd) {

  batch_list <- as.list(as.character(unique(batch)))
  batch_oneType <- unlist(batch_list)[which(kmeansK == 1)]
  batch_num <- table(batch)[as.character(unique(batch))]

  # if(ncol(exprs_mat)>=5000){ rpca_q = 0 }else if(ncol(exprs_mat)>=2000){ rpca_q =1 }else{ rpca_q =2 }

  #############################

  # if (parallel) {
    pca <- BiocParallel::bplapply(batch_list, function(this_batch_list) {
      computePCA_byHVGMarker(this_batch_list = this_batch_list, batch = batch, batch_oneType = batch_oneType,
                             marker = marker, exprs_mat = exprs_mat, HVG_list = HVG_list, fast_svd = fast_svd)
    }, BPPARAM = parallelParam)
  # } else {
  #   pca <- lapply(batch_list, function(this_batch_list) {
  #     computePCA_byHVGMarker(this_batch_list = this_batch_list, batch = batch, batch_oneType = batch_oneType,
  #                            marker = marker, exprs_mat = exprs_mat, HVG_list = HVG_list, fast_svd = fast_svd)
  #   })
  # }

  #############################

  names(pca) <- unlist(batch_list)

  clustering_res <- list()
  clustering_res_pt_dist <- list()

  # if (parallel) {
    res <- BiocParallel::bplapply(seq_len(length(pca)), function(j) {
      pca_current <- pca[[j]]
      if (!is.null(pca_current)) {
        kmeans_res <- stats::kmeans(pca_current[, 1:10], centers = kmeansK[j], nstart = 1000)
        clustering_res_tmp <- kmeans_res$cluster
        # exprs_current <- exprs_mat[HVG_list[[j]], batch == batch_list[j]]
        if (!is.null(marker)) {
          exprs_current <- exprs_mat[marker, batch == batch_list[j]]
        } else {
          exprs_current <- exprs_mat[HVG_list[[j]], batch == batch_list[j]]
        }
        clustering_res_pt_dist_tmp <- lapply(1:kmeansK[j], function(y) {

          centroidDist(exprs_current[, kmeans_res$cluster == y, drop = FALSE])
        })
        clustering_res_pt_dist_tmp <- unlist(clustering_res_pt_dist_tmp)
        clustering_res_pt_dist_tmp <- clustering_res_pt_dist_tmp[names(clustering_res_tmp)]
      } else {
        clustering_res_tmp <- rep(1, batch_num[j])
        names(clustering_res_tmp) <- colnames(exprs_mat[, batch == batch_list[j]])
        if (!is.null(marker)) {
          clustering_res_pt_dist_tmp <- centroidDist(exprs_mat[marker, batch == batch_list[j]])
        } else {
          clustering_res_pt_dist_tmp <- centroidDist(exprs_mat[HVG_list[[j]], batch == batch_list[j]])
        }
      }
      list(clustering_res_tmp, clustering_res_pt_dist_tmp)
    }, BPPARAM = parallelParam)

    clustering_res <- lapply(res, function(x) x[[1]])
    clustering_res_pt_dist <- lapply(res, function(x) x[[2]])
  # } else {
  #   for (j in 1:length(pca)) {
  #     pca_current <- pca[[j]]
  #     if (!is.null(pca_current)) {
  #       kmeans_res <- stats::kmeans(pca_current[, 1:10], centers = kmeansK[j], nstart = 1000)
  #       clustering_res[[j]] <- kmeans_res$cluster
  #       # exprs_current <- exprs_mat[HVG_list[[j]], batch == batch_list[j]]
  #       if (!is.null(marker)) {
  #         exprs_current <- exprs_mat[marker, batch == batch_list[j]]
  #       } else {
  #         exprs_current <- exprs_mat[HVG_list[[j]], batch == batch_list[j]]
  #       }
  #       clustering_res_pt_dist[[j]] <- lapply(1:kmeansK[j], function(y) {
  #
  #         centroidDist(exprs_current[, kmeans_res$cluster == y, drop = FALSE])
  #       })
  #       clustering_res_pt_dist[[j]] <- unlist(clustering_res_pt_dist[[j]])
  #       clustering_res_pt_dist[[j]] <- clustering_res_pt_dist[[j]][names(clustering_res[[j]])]
  #     } else {
  #       clustering_res[[j]] <- rep(1, batch_num[j])
  #       names(clustering_res[[j]]) <- colnames(exprs_mat[, batch == batch_list[j]])
  #       if (!is.null(marker)) {
  #         clustering_res_pt_dist[[j]] <- centroidDist(exprs_mat[marker, batch == batch_list[j]])
  #       } else {
  #         clustering_res_pt_dist[[j]] <- centroidDist(exprs_mat[HVG_list[[j]], batch == batch_list[j]])
  #       }
  #     }
  #   }
  # }


  return(list(clustering_list = clustering_res, clustering_distProp = clustering_res_pt_dist))
}
######################################################################################################
centroidDist <- function(exprs_mat) {
  centroid_batch <- matrixStats::rowMedians(exprs_mat)
  point_dist <- colSums((exprs_mat - centroid_batch)^2)
  point_rank <- rank(point_dist)
  point_rank <- point_rank/length(point_rank)
  return(point_rank)
}

###################################################################################################### Function to find the mutual nearest clusters
compute_dist_mat_med <- function(k, exprs_mat, clustering_list, combine_pair, dist) {
  ## We go through every pairwise batches print(k) Extract the cell type information and the expression
  ## matrices
  res1 <- clustering_list[[combine_pair[1, k]]]
  res2 <- clustering_list[[combine_pair[2, k]]]
  mat1 <- exprs_mat[, names(which(res1 == 1))]
  mat2 <- exprs_mat[, names(which(res2 == 1))]

  ## The distance between matrices are calculated as such...
  if (dist == "cosine") {
    dist_mat <- dist(t(mat1), t(mat2), method = "cosine")
  } else if (dist == "cor") {
    dist_mat <- 1 - cor((mat1), (mat2))
  } else {
    dist_mat <- pdist::pdist(t(mat1), t(mat2))
  }


  dist_mat <- as.matrix(dist_mat)

  ## The dist_res (distance measure between batches) is then the median of all pairwise distances
  return(median(dist_mat))
}

compute_dist_res <- function(i, res1, res2, exprs_mat, dist, dist_res) {
  res_tmp <- c()
  for (j in 1:max(res2)) {
    mat1 <- exprs_mat[, names(which(res1 == i))]
    mat2 <- exprs_mat[, names(which(res2 == j))]
    if (dist == "cosine") {
      dist_mat <- dist(t(mat1), t(mat2), method = "cosine")
    } else if (dist == "cor") {
      dist_mat <- 1 - cor((mat1), (mat2))
    } else {
      dist_mat <- pdist::pdist(t(mat1), t(mat2))
    }
    dist_mat <- as.matrix(dist_mat)
    res_tmp <- c(res_tmp, median(dist_mat))
  }
  dist_res[[i]] <- res_tmp
}

findMNC <- function(exprs_mat, clustering_list, dist = "euclidean", parallelParam) {

  batch_num <- length(clustering_list)
  names(clustering_list) <- paste("Batch", 1:batch_num, sep = "")

  ## Check which batch has only one cluster
  batch_oneType <- which(unlist(lapply(clustering_list, function(x) length(levels(as.factor(x))) ==
                                         1)))
  ## If there existsome batch_oneType
  if (length(batch_oneType) != 0) {
    ## And if all batch_oneType == num of batches, i.e. every batch only contains one cell type
    if (length(batch_oneType) == batch_num) {
      combine_pair <- utils::combn(batch_num, 2)
      batch_oneType <- NULL
      allones <- TRUE
    } else {
      ## if at least some batch contains more than 1 cell type Then take away the batches with only one cell
      ## type and then iterate through all combn
      combine_pair <- combn(c(1:batch_num)[-batch_oneType], 2)

      ## And then for those batches with only one cell type, we bind to the previous generated combn
      for (i in batch_oneType) {
        for (j in c(1:batch_num)[-batch_oneType]) {
          combine_pair <- cbind(combine_pair, c(i, j))
        }
      }
      allones <- FALSE
    }
  } else {
    combine_pair <- utils::combn(batch_num, 2)
    allones <- FALSE
  }




  mnc <- list()

  ## If there are only two batches containing only two cell types, then finding MNN is trivial. Return
  ## NULL
  if (allones & batch_num == 2) {
    return(NULL)
  } else if (allones) {
    ## If every batch contains only one cell type...
    dist_res <- matrix(NA, nrow = batch_num, ncol = batch_num)



    # if (parallel) {
      dist_mat_med <- BiocParallel::bplapply(seq_len(ncol(combine_pair)), function(k) {
        compute_dist_mat_med(k = k, exprs_mat = exprs_mat, clustering_list = clustering_list,
                             combine_pair = combine_pair, dist = dist)
      }, BPPARAM = parallelParam)
    # } else {
    #   dist_mat_med <- lapply(seq_len(ncol(combine_pair)), function(k) {
    #     compute_dist_mat_med(k = k, exprs_mat = exprs_mat, clustering_list = clustering_list,
    #                          combine_pair = combine_pair, dist = dist)
    #   })
    # }


    for (k in 1:ncol(combine_pair)) {
      dist_res[combine_pair[1, k], combine_pair[2, k]] <- dist_res[combine_pair[2, k], combine_pair[1,
                                                                                                    k]] <- dist_mat_med[[k]]
    }

    ## The neighbour_res is then which ever two pairs of batches that are closes to each other
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
      colnames(mnc[[i]]) <- c(paste("Batch", mnc_mat[i, 1], sep = ""), paste("Batch", mnc_mat[i,
                                                                                              2], sep = ""))
    }
  } else {
    for (k in 1:ncol(combine_pair)) {
      dist_res <- list()
      # print(k)
      res1 <- clustering_list[[combine_pair[1, k]]]
      res2 <- clustering_list[[combine_pair[2, k]]]


      # if (parallel) {
        dist_res <- BiocParallel::bplapply(seq_len(max(res1)), function(i) {
          compute_dist_res(i = i, res1 = res1, res2 = res2, exprs_mat = exprs_mat, dist = dist,
                           dist_res)
        }, BPPARAM = parallelParam)
      # } else {
      #   dist_res <- lapply(seq_len(max(res1)), function(i) {
      #     compute_dist_res(i = i, res1 = res1, res2 = res2, exprs_mat = exprs_mat, dist = dist,
      #                      dist_res)
      #   })
      # }


      dist_res <- do.call(rbind, dist_res)
      neighbour_batch1 <- apply(dist_res, 1, which.min)
      neighbour_batch2 <- apply(dist_res, 2, which.min)
      mnc_tmp <- c()
      for (l in seq_len(length(neighbour_batch1))) {
        if (neighbour_batch2[neighbour_batch1[l]] == l) {
          mnc_tmp <- rbind(mnc_tmp, c(l, neighbour_batch1[l]))
        }
      }
      mnc[[k]] <- mnc_tmp
      colnames(mnc[[k]]) <- c(paste("Batch", combine_pair[1, k], sep = ""), paste("Batch", combine_pair[2,
                                                                                                        k], sep = ""))
    }
  }  ## End else

  ############################################################### Function to perform network analysis


  edge_list <- do.call(rbind, lapply(mnc, function(x) t(apply(x, 1, function(y) paste(colnames(x),
                                                                                      y, sep = "_")))))

  if (is.null(edge_list)) {
    return(NULL)
  } else {
    g <- igraph::graph_from_edgelist(edge_list, directed = FALSE)
    igraph::plot.igraph(g)
    mnc <- igraph::fastgreedy.community(g)
    mnc_df <- data.frame(group = as.numeric(mnc$membership), batch = as.numeric(gsub("Batch", "",
                                                                                     gsub("_.*", "", mnc$names))), cluster = as.numeric(gsub(".*_", "", mnc$names)))

    if (allones) {
      mnc_df_new <- mnc_df
      batch_oneType <- c(1:batch_num)[-c(mnc_mat)]
      for (i in batch_oneType) {
        print(i)
        neighbour_order <- rank(dist_res[i, ], na.last = TRUE)
        group_order1 <- mnc_df[mnc_df[, "batch"] == which(neighbour_order == 1), "group"]
        group_order2 <- mnc_df[mnc_df[, "batch"] == which(neighbour_order == 2), "group"]
        if (group_order1 == group_order2) {
          mnc_df_new <- rbind(mnc_df_new, c(group_order1, i, 1))
        }
      }
      mnc_df <- mnc_df_new
    }

    return(mnc_df)
  }

}



###################################################################################################### Function to create replicate based on the mutual nearest cluster results
mncReplicate <- function(clustering_list, clustering_distProp, replicate_prop, mnc_df) {
  batch_num <- length(clustering_list)
  if (!is.null(mnc_df)) {
    idx_noRep <- list()

    # For each batches, check whether there is any clusters that do not have replicates for(j in
    # 1:length(clustering_list)){ idx_noRep[[j]]<-which(!(1:max(clustering_list[[j]]))%in%mnc[,j]) }
    for (j in 1:length(clustering_list)) {
      idx_noRep[[j]] <- which(!(1:max(clustering_list[[j]])) %in% mnc_df[mnc_df$batch == j, "cluster"])
    }

    replicate_vector <- rep(NA, length(unlist(clustering_list)))
    names(replicate_vector) <- names(unlist(clustering_list))
    clustering_distProp <- unlist(clustering_distProp)
    replicate_size <- table(mnc_df$group)
    for (i in 1:max(mnc_df$group)) {
      tmp_names <- c()

      # For each batch in this replicate
      mnc_df_sub <- mnc_df[mnc_df$group == i, ]
      for (l in 1:replicate_size[i]) {
        # tmp_names<-c(tmp_names,names(which(clustering_list[[l]]==mnc[i,l])))

        tmp_names <- c(tmp_names, names(which(clustering_list[[mnc_df_sub[l, "batch"]]] == mnc_df_sub[l,
                                                                                                      "cluster"])))
      }

      replicate_vector[tmp_names[clustering_distProp[tmp_names] <= replicate_prop]] <- paste("Replicate",
                                                                                             i, sep = "_")
    }

    current_idx <- max(mnc_df$group)
    for (j in 1:length(clustering_list)) {
      if (length(idx_noRep[[j]]) != 0) {
        for (k in 1:length(idx_noRep[[j]])) {
          tmp_names <- names(which(clustering_list[[j]] == idx_noRep[[j]][k]))
          replicate_vector[tmp_names[clustering_distProp[tmp_names] <= replicate_prop]] <- paste("Replicate",
                                                                                                 current_idx + k, sep = "_")
        }
        current_idx <- current_idx + length(idx_noRep[[j]])
      }
    }
    replicate_vector[is.na(replicate_vector)] <- 1:sum(is.na(replicate_vector))
  } else {
    replicate_vector <- rep(NA, length(unlist(clustering_list)))
    names(replicate_vector) <- names(unlist(clustering_list))
    clustering_distProp <- unlist(clustering_distProp)
    current_idx <- 1
    for (j in 1:length(clustering_list)) {
      for (k in 1:max(clustering_list[[j]])) {
        tmp_names <- names(which(clustering_list[[j]] == k))
        replicate_vector[tmp_names[clustering_distProp[tmp_names] <= replicate_prop]] <- paste("Replicate",
                                                                                               current_idx + k, sep = "_")
      }
      current_idx <- current_idx + max(clustering_list[[j]])
    }
    replicate_vector[is.na(replicate_vector)] <- 1:sum(is.na(replicate_vector))
  }

  return(replicate_vector)
}

###################################################################################################### Function to create replicates based on wanted varition
wvReplicate <- function(exprs_mat, WV, WV_marker, replicate_vector) {
  names(WV) <- colnames(exprs_mat)
  if (!is.null(WV_marker)) {
    marker_expr <- lapply(WV_marker, function(x) aggregate(exprs_mat[x, ], by = list(replicate_vector),
                                                           FUN = mean))
    names(marker_expr) <- WV_marker

    marker_expr <- do.call(cbind, lapply(marker_expr, function(x) x[, 2]))
    rownames(marker_expr) <- names(table(replicate_vector))
    km <- stats::kmeans(marker_expr, centers = 2)
    tab_max <- table(apply(km$centers, 2, which.max))
    marker_cluster <- names(tab_max)[which.max(tab_max)]
    marker_replicate <- names(which(km$cluster == marker_cluster))
    replicate_vector[replicate_vector %in% marker_replicate] <- paste(replicate_vector[replicate_vector %in%
                                                                                         marker_replicate], WV[replicate_vector %in% marker_replicate], sep = "_")

    # repMat_new <- ruv::replicate.matrix(as.factor(replicate_vector))
  } else {
    replicate_vector <- paste(replicate_vector, WV, sep = "_")
    # repMat_new <- ruv::replicate.matrix(as.factor(paste(replicate_vector, WV, sep = '_')))
  }

  return(replicate_vector)
}

###################################################################################################### Function to create replicates based on known cell types
supervisedReplicate <- function(exprs_mat, cell_type, replicate_prop) {
  clustering_distProp <- lapply(unique(cell_type), function(x) centroidDist(exprs_mat[, cell_type ==
                                                                                        x]))
  clustering_distProp <- unlist(clustering_distProp)[names(cell_type)]
  replicate_vector <- rep(NA, length(cell_type))
  replicate_vector[clustering_distProp <= replicate_prop] <- cell_type[clustering_distProp <= replicate_prop]
  replicate_vector[is.na(replicate_vector)] <- 1:sum(is.na(replicate_vector))
  return(replicate_vector)
}


###################################################################################################### Function to create replicates based on known cell types
computePCA_byHVGMarker <- function(this_batch_list, batch, batch_oneType, marker, exprs_mat, HVG_list,
                                   fast_svd) {
  ## If there is kmeansK == 1 in a batch, then we will return NULL
  if (this_batch_list %in% batch_oneType) {
    return(NULL)
  } else {
    ## If there is kmeansK > 1 in a batch, then we will return compute PCAs

    ## If there exist some user-supplied markers, then we will use them to perform PCA.  Otherwise we use
    ## the HVGs computed in each batch
    if (!is.null(marker)) {
      sub_exprs_mat <- t(exprs_mat[marker, batch == this_batch_list])
    } else {
      sub_exprs_mat <- t(exprs_mat[HVG_list[[this_batch_list]], batch == this_batch_list])
    }

    ## If fast_svd option was enabled, then we will use irlba to speed up the calculations.
    if (fast_svd) {
      result <- irlba::prcomp_irlba(x = sub_exprs_mat, n = 10, scale. = TRUE, maxit = 1000)$x
    } else {
      result <- prcomp(sub_exprs_mat, scale. = TRUE)$x
    }
  }
  rownames(result) <- rownames(sub_exprs_mat)
  return(result)
}
