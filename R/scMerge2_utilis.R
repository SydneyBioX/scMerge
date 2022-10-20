#' @importFrom scran getTopHVGs modelGeneVar

feature_selection <- function(exprsMat, batch, top_n = 2000, use_bpparam = SerialParam()) {
    combined.dec <- scran::modelGeneVar(exprsMat, block = batch, BPPARAM = use_bpparam)
    chosen.hvgs <- scran::getTopHVGs(combined.dec, n = top_n)
    return(chosen.hvgs)
}





#' @import BiocParallel
#' @importFrom BiocSingular runSVD
#' @importFrom DelayedArray DelayedArray t sweep colMeans

pseudoRUVIII <- function(Y, Y_pseudo, M, ctl, k = NULL, eta = NULL, 
                         include.intercept = TRUE, inputcheck = TRUE, 
                         fullalpha = NULL, return.info = FALSE, 
                         subset = NULL,
                         BPPARAM = BiocParallel::SerialParam(), 
                         BSPARAM = BiocSingular::ExactParam(),
                         normalised = TRUE) {
    
    
    
    
    m = nrow(Y)
    n = ncol(Y)
    
    if (is.data.frame(Y_pseudo)) {
        Y = data.matrix(Y_pseudo)
    }
    m_pseudo = nrow(Y_pseudo)
    n_pseudo = ncol(Y_pseudo)
    
    
    M = replicate.matrix(M)
    ctl = tological(ctl, n)
    if (inputcheck) {
        # if (m > n) 
        #   warning("m is greater than n!  This is not a problem itself,
        #           but may indicate that you need to transpose your data matrix.  
        #           Please ensure that rows correspond to observations (e.g. microarrays) 
        #           and columns correspond to features (e.g. genes).")
        if (sum(is.na(Y)) > 0) 
            warning("Y contains missing values.  This is not supported.")
        if (sum(Y == Inf, na.rm = TRUE) + sum(Y == -Inf, na.rm = TRUE) > 0) 
            warning("Y contains infinities.  This is not supported.")
    }
    
    Y_pseudo = RUV1(Y_pseudo, eta, ctl, include.intercept = include.intercept)
    
    
    
    
    if (class(BSPARAM) != "ExactParam") {
        svd_k <- k
        svd_k <- min(m_pseudo - ncol(M), sum(ctl), svd_k, na.rm = TRUE)
    } else {
        svd_k <- min(m_pseudo - ncol(M), sum(ctl), na.rm = TRUE)
    }
    
    ## m represent the number of samples/observations ncol(M)
    ## represent the number of replicates If the replicate matrix
    ## is such that we have more replicates than samples, then
    ## RUV3 is not appropriate, thus, we return the Original input
    ## matrix
    if (ncol(M) >= m_pseudo | k == 0) {
        newY <- Y_pseudo
        fullalpha <- NULL
        return(newY)
    } 
    
    if (is.null(fullalpha)) {
        ## The main RUVIII process Applies the residual operator of a
        ## matrix M to a matrix Y Y0 has the same dimensions as Y,
        ## i.e. m rows (observations) and n columns (genes).
        
        Y0 <- my_residop(Y_pseudo, M)
        svdObj <- BiocSingular::runSVD(x = Y0, k = svd_k, BPPARAM = BPPARAM, BSPARAM = BSPARAM)
        
        fullalpha <- DelayedArray::t(svdObj$u[, seq_len(svd_k), drop = FALSE]) %*% Y_pseudo
    }  ## End is.null(fullalpha)
    
    ###############
    alpha <- fullalpha[seq_len(min(k, nrow(fullalpha))), , drop = FALSE]
    
    alpha <- DelayedArray::DelayedArray(alpha)
    # newY <-  Y - W %*% alpha
    
    if (normalised) {
        
        ac <- alpha[, ctl, drop = FALSE]
        
        Y <- DelayedArray::DelayedArray(Y)
        Y_stand <- DelayedArray::sweep(Y, 2, DelayedArray::colMeans(Y), "-")
        print(dim(Y_stand))
        W <- Y_stand[, ctl] %*% DelayedArray::t(ac) %*% solve(ac %*% DelayedArray::t(ac))
        W <- DelayedArray::DelayedArray(W)
        
        
        if (!is.null(subset)) {
            newY <- Y[, subset] - W %*% alpha[, subset]
        } else {
            newY <- Y - W %*% alpha
        }
        
        if (!return.info) {
            return(newY)
        } else {
            return(list(newY = newY, M = M, fullalpha = fullalpha))
        }
    } else {
        return(list(M = M, fullalpha = fullalpha))
    }
    
}



#' @importFrom scran buildSNNGraph
#' @importFrom igraph cluster_louvain
#' @importFrom BiocParallel bplapply
#' @importFrom batchelor cosineNorm

group_wtihin_batch <- function(exprsMat = NULL,
                               cellTypes = NULL,
                               pseudobulk_batch_list = NULL,
                               pseudobulk_batch = NULL,
                               pseudobulk_sample_list = NULL,
                               pseudobulk_sample = NULL,
                               chosen.hvg = rownames(exprsMat),
                               k_celltype = 10,
                               use_bnparam = BiocNeighbors::AnnoyParam(),
                               use_bpparam = BiocParallel::SerialParam(),
                               seed = 1,
                               verbose = TRUE) {
    
    if (is.null(cellTypes)) {
        
        if (verbose) {
            print("Cluster within batch")
        }
        clust <- BiocParallel::bplapply(pseudobulk_batch_list, function(x) {
            set.seed(seed)
            g <- scran::buildSNNGraph(exprsMat[chosen.hvg, pseudobulk_batch == x], 
                                      k = k_celltype, 
                                      BNPARAM = use_bnparam)
            res <- igraph::cluster_louvain(g)$membership
            names(res) <- colnames(exprsMat)[pseudobulk_batch == x]
            res
        }, BPPARAM = use_bpparam)
        clust <- unlist(clust)
        clust <- clust[colnames(exprsMat)]
        clust <- split(clust, pseudobulk_sample)
        clust <- lapply(clust, function(x) {
            as.numeric(factor(x))
        })
        
    } else {
        
        if (verbose) {
            print("Use existing cell type info within batch")
        }
        
        clust <- BiocParallel::bplapply(pseudobulk_sample_list, function(x) {
            res <- as.numeric(factor(cellTypes[pseudobulk_sample == x]))
            names(res) <- colnames(exprsMat)[pseudobulk_sample == x]
            res
        }, BPPARAM = use_bpparam)
    }
    
    return(clust)
}

#' @importFrom methods as

construct_pseudoBulk <- function(exprsMat = NULL,
                                 exprsMat_counts = NULL,
                                 pseudobulk_sample = NULL,
                                 pseudobulk_sample_list = NULL,
                                 clust = NULL,
                                 cosineNorm = NULL,
                                 pseudoBulk_fn = "create_pseudoBulk",
                                 k_psuedoBulk = 30, 
                                 use_bpparam = BiocParallel::SerialParam(),
                                 verbose = TRUE) {
    
    aggExprs <- list()
    for(i in seq_along(pseudobulk_sample_list)){
        if (pseudoBulk_fn == "create_pseudoBulk_divide") {
            res <- create_pseudoBulk_divide(exprsMat_counts[, pseudobulk_sample == pseudobulk_sample_list[i]],
                                            clust[[i]], 
                                            k_fold = k_psuedoBulk,
                                            use_bpparam = use_bpparam)
        }
        
        
        if (pseudoBulk_fn == "create_pseudoBulk_pool_divide") {
            res <- create_pseudoBulk_pool_divide(exprsMat_counts[, pseudobulk_sample == pseudobulk_sample_list[i]],
                                                 clust[[i]], 
                                                 k_fold = k_psuedoBulk,
                                                 use_bpparam = use_bpparam)
        }
        
        
        if (pseudoBulk_fn == "create_pseudoBulk") {
            res <- create_pseudoBulk(exprsMat[, pseudobulk_sample == pseudobulk_sample_list[i]],
                                     clust[[i]], 
                                     k_fold = k_psuedoBulk,
                                     use_bpparam = use_bpparam)
        }
        
        rownames(res) <- paste(i, rownames(res), sep = "_")
        aggExprs[[i]] <-res
    }
    
    
    
    bulkExprs <- t(do.call(rbind, aggExprs))
    
    
    if (as.character(substitute(pseudoBulk_fn)) %in% c("create_pseudoBulk_pool_divide",
                                                       "create_pseudoBulk_divide") & cosineNorm) {
        bulkExprs <- batchelor::cosineNorm(bulkExprs, BPPARAM = use_bpparam)
    }
    
    
    rownames(bulkExprs) <- rownames(exprsMat)
    
    bulkExprs[is.infinite(bulkExprs)] <- 0
    bulkExprs <- methods::as(bulkExprs, "CsparseMatrix")
    
    pseudo_batch <- unlist(lapply(strsplit(colnames(bulkExprs), "_"), "[[", 1))
    
    
    if (verbose) {
        cat("Dimension of psuedo-bulk expression: ")
        print(dim(bulkExprs))
    }
    
    pseudo_batch <- unlist(lapply(strsplit(colnames(bulkExprs), "_"), "[[", 1))
    
    
    
    bulk_clustering_res <- lapply(aggExprs, function(x) {
        res <- unlist(lapply(strsplit(rownames(x), "_"), "[[", 2))
        res <- as.numeric(res)
        names(res) <- rownames(x)
        res
    })
    
    bulk_clustering_distProp <- lapply(bulk_clustering_res, function(x) {
        res <- rep(1, length(x))
        names(res) <- names(x)
        res
    })
    
    
    return(list(bulkExprs = bulkExprs,
                pseudo_batch = pseudo_batch,
                bulk_clustering_res = bulk_clustering_res,
                bulk_clustering_distProp = bulk_clustering_distProp))
}





#' @importFrom ruv replicate.matrix
#' @importFrom methods as is

aggregate.Matrix <- function(x, groupings=NULL) {
    if (!methods::is(x,'Matrix')) {
        x <- methods::as(as.matrix(x), "CsparseMatrix")
    }
    
    groupings2 <- paste("A", groupings, sep = "")
    
    if (length(unique(groupings2)) > 1) {
        
        mapping <- methods::as(ruv::replicate.matrix(groupings2), "CsparseMatrix")
        colnames(mapping) <- substring(colnames(mapping), 2)
        mapping <- mapping[, levels(factor(groupings))]
        
    } else {
        mapping <- methods::as(matrix(rep(1, length(groupings2)), ncol = 1), "CsparseMatrix")
        colnames(mapping) <- unique(groupings)
    }
    
    result <- t(mapping) %*% x
    
    return(result)
}




#' @importFrom scater normalizeCounts
#' @importFrom BiocParallel bplapply


create_pseudoBulk_pool_divide <- function(exprsMat_counts, 
                                          cell_info, 
                                          k_fold = 30,
                                          use_bpparam = BiocParallel::SerialParam()) {
    
    
    cell_info <- as.character(cell_info)
    exprsMat_pseudo <- BiocParallel::bplapply(split(seq_len(ncol(exprsMat_counts)), cell_info), function(x) {
        if (length(x) > k_fold) {
            # pool
            total_counts <- colSums(exprsMat_counts[, x])
            bins <- cut(rank(total_counts), k_fold)
            res <- aggregate.Matrix(t(exprsMat_counts[, x, drop = FALSE]), 
                                    bins)
            bins_tab <- table(bins)
            cellTypes_n_mat <- matrix(rep(bins_tab[rownames(res)], 
                                          ncol(t(exprsMat_counts[, x]))), 
                                      nrow = nrow(res), byrow = FALSE)
            res <- lapply(seq_len(nrow(res)), function(k) {
                stats::rbinom(ncol(res), round(res[k, ]), 
                              1/cellTypes_n_mat[k, ])
            })
            res <- do.call(cbind, res)
            res <- scater::normalizeCounts((res))
            colnames(res) <- seq_len(ncol(res))
            res
        } else {
            res <- exprsMat_counts[, x, drop = FALSE]
            res <- scater::normalizeCounts((res), log = TRUE)
            colnames(res) <- seq_len(ncol(res))
            res
        }
        
        
    }, BPPARAM = use_bpparam)
    
    exprsMat_pseudo <- lapply(1:length(exprsMat_pseudo), function(i) {
        colnames(exprsMat_pseudo[[i]]) <- paste(i, colnames(exprsMat_pseudo[[i]]), sep = "_")
        exprsMat_pseudo[[i]]
    })
    
    
    exprsMat_pseudo <- do.call(cbind, exprsMat_pseudo)
    
    return(t(exprsMat_pseudo))
}


#' @importFrom scater normalizeCounts
#' @importFrom BiocParallel bplapply
#' @importFrom cvTools cvFolds

create_pseudoBulk_divide <- function(exprsMat_counts, 
                                     cell_info, 
                                     k_fold = 30, 
                                     use_bpparam = BiocParallel::SerialParam()) {
    
    
    
    
    exprsMat_pseudo <- BiocParallel::bplapply(split(seq_len(ncol(exprsMat_counts)), cell_info), function(x) 
    {
        if (length(x) > k_fold) {
            cv <- cvTools::cvFolds(length(x), K = k_fold)
            bins <- cv$which[cv$subsets]
            res <- aggregate.Matrix(t(exprsMat_counts[, x, drop = FALSE]), 
                                    bins)
            bins_tab <- table(bins)
            cellTypes_n_mat <- matrix(rep(bins_tab, 
                                          ncol(t(exprsMat_counts[, x]))), 
                                      nrow = length(bins_tab), byrow = FALSE)
            res <- lapply(seq_len(nrow(res)), function(k) stats::rbinom(ncol(res),  round(res[k, ]), 
                                                                        1/cellTypes_n_mat[k, ]))
            res <- do.call(cbind, res)
            res <- scater::normalizeCounts((res))
            colnames(res) <- as.numeric(factor(names(bins_tab)))
            res
        } else {
            res <- exprsMat_counts[, x, drop = FALSE]
            colnames(res) <- seq_len(ncol(res))
            res
        }
        
        
    }, BPPARAM = use_bpparam)
    
    
    exprsMat_pseudo <- lapply(1:length(exprsMat_pseudo), function(i) {
        colnames(exprsMat_pseudo[[i]]) <- paste(i, colnames(exprsMat_pseudo[[i]]), sep = "_")
        exprsMat_pseudo[[i]]
    })
    
    
    exprsMat_pseudo <- do.call(cbind, exprsMat_pseudo)
    
    
    return(t(exprsMat_pseudo))
}



#' @importFrom scater normalizeCounts
#' @importFrom BiocParallel bplapply
#' @importFrom cvTools cvFolds


create_pseudoBulk <- function(exprsMat, cell_info, k_fold = 30, use_bpparam = BiocParallel::SerialParam()) {
    
    k_fold <- min(ncol(exprsMat), k_fold)
    cv <- cvTools::cvFolds(ncol(exprsMat), K = k_fold)
    
    exprsMat_pseudo <- BiocParallel::bplapply(seq_len(k_fold), function(i) {
        subset_idx <- cv$subsets[cv$which == i]
        cellType_tab <- table(droplevels(factor(cell_info[subset_idx])))
        cellTypes_n_mat <- matrix(rep(cellType_tab, 
                                      nrow(exprsMat)), 
                                  nrow = length(cellType_tab), byrow = FALSE)
        rownames(cellTypes_n_mat) <- names(cellType_tab)
        
        res <- aggregate.Matrix(t(exprsMat[, subset_idx]), 
                                cell_info[subset_idx])
        cellTypes_n_mat <- cellTypes_n_mat[rownames(res), ]
        res <- res/cellTypes_n_mat
        
        rownames(res) <- paste(rownames(res), i, sep = "_")
        res
    }, BPPARAM = use_bpparam)
    
    exprsMat_pseudo <- do.call(rbind, exprsMat_pseudo)
    
    return(exprsMat_pseudo)
}






check_input2 <- function(exprsMat, batch, 
                         cellTypes, condition, 
                         ctl, chosen.hvg, return_subset_genes,
                         exprsMat_counts){
    
    #### Checking input exprsMat
    
    if (is.null(exprsMat)) {
        stop("The 'exprsMat' argument is NULL.")
    }
    
    
    #### Checking if the cell names are non-unique
    cell_names <- colnames(exprsMat)
    
    if (length(cell_names) != length(unique(cell_names)) | is.null(cell_names)) {
        stop("Column names of the input exprsMat object must not contain duplicates nor NULL")
    }
    
    
    
    
    #### Check batch input
    
    if (is.null(batch)) {
        stop("The 'batch' argument is NULL.")
    }
    
    if (any(is.na(batch))) {
        stop("NA's found the batch column, please remove")
    }
    
    
    if (ncol(exprsMat) != length(batch)) {
        stop("The length of batch is not equal to the number of column in exprsMat.")
    }
    
    #### Check cell types input
    
    if (!is.null(cellTypes)) {
        
        if (any(is.na(cellTypes))) {
            stop("NA's found the cellTypes column, please remove")
        }
        
        
        if (ncol(exprsMat) != length(cellTypes)) {
            stop("The length of cellTypes is not equal to the number of column in exprsMat.")
        }
    }
    

    
    #### Check condition input
    
    if (!is.null(condition)) {
        
        if (any(is.na(condition))) {
            stop("NA's found the condition column, please remove.")
        }
        
        
        if (ncol(exprsMat) != length(condition)) {
            stop("The length of condition is not equal to the number of column in exprsMat.")
        }
    }
    
    
    #### check ctl, chosen.hvg and return_subset_genes
    
    if (sum(ctl %in% rownames(exprsMat)) == 0) {
        stop("There is no ctl genes found in the rownames of exprsMat.")
    }
    


    
    if (!is.null(chosen.hvg)) {
        if (sum(chosen.hvg %in% rownames(exprsMat)) == 0) {
            stop("There is no chosen.hvg genes found in the rownames of exprsMat.")
        }
    }
    
    if (!is.null(return_subset_genes)) {
        if (sum(return_subset_genes %in% rownames(exprsMat)) == 0) {
            stop("There is no return_subset_genes genes found in the rownames of exprsMat.")
        }
    }

    
    #### Check cell types input
    
    if (!is.null(exprsMat_counts)) {
        
        if (ncol(exprsMat_counts) != ncol(exprsMat)) {
            stop("The number of column in exprsMat_counts is not equal to the number of column in exprsMat")
        }
        
    }
    

    

}

