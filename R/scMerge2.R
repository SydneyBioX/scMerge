#' @title Perform the scMerge2 algorithm 
#' @description Merge single-cell data from different batches and experiments 
#' leveraging (pseudo)-replicates, control genes and pseudo-bulk.
#' @author Yingxin Lin
#' @param exprsMat A gene (row) by cell (column) log-transformed matrix to be adjusted.
#' @param batch A vector indicating the batch information for each cell in the batch-combined matrix.
#' @param cellTypes An optional vector indicating the cell type information for each cell in the batch-combined matrix. 
#' If it is \code{NULL}, pseudo-replicate procedure will be run to identify cell type.
#' @param condition An optional vector indicating the condition information for each cell in the batch-combined matrix. 
#' @param ctl A character vector of negative control. It should have a non-empty intersection with the rows of exprsMat
#' @param chosen.hvg An optional character vector of highly variable genes chosen.
#' @param ruvK An integer indicates the number of unwanted variation factors that are removed, default is 20.
#' @param use_bpparam A \code{BiocParallelParam} class object from the \code{BiocParallel} package is used. Default is SerialParam().
#' @param use_bsparam A \code{BiocSingularParam} class object from the \code{BiocSingular} package is used. Default is RandomParam().
#' @param use_bnparam A \code{BiocNeighborsParam} class object from the \code{BiocNeighbors} package is used. Default is AnnoyParam().
#' @param pseudoBulk_fn A character indicates the way of pseudobulk constructed.
#' @param k_pseudoBulk An integer indicates the number of pseudobulk constructed within each cell grouping. Default is 30.
#' @param k_celltype An integer indicates the number of nearest neighbours used in \code{buildSNNGraph} when grouping cells within each batch. Default is 10.
#' @param exprsMat_counts A gene (row) by cell (column) counts matrix to be adjusted.
#' @param cosineNorm A logical vector indicates whether cosine normalisation is performed on input data.
#' @param return_subset If \code{TRUE}, adjusted matrix of only a subset of genes (hvg or indicates in \code{return_subset_genes}) will be return. 
#' @param return_subset_genes An optional character vector of indicates the subset of genes will be adjusted.
#' @param return_matrix A logical vector indicates whether the adjusted matrix is calculated and returned. 
#' If \code{FALSE}, then only the estimated parameters will be returned.
#' @param verbose If \code{TRUE}, then all intermediate steps will be shown. Default to \code{FALSE}.
#' @param seed A numeric input indicates the seed used.
#' 
#' @import SummarizedExperiment
#' @importFrom S4Vectors metadata
#' @importFrom BiocParallel SerialParam bpparam
#' @importFrom BiocSingular ExactParam
#' @importFrom BiocNeighbors AnnoyParam
#' @importFrom batchelor cosineNorm
#' @importFrom DelayedMatrixStats rowSums2 colSums2
#' @importFrom methods as is
#' 
#' @return Returns a \code{list} object with following components:
#' \itemize{
#' \item{newY: }{if \code{return_matrix} is \code{TRUE}, the adjusted matrix will be return.}
#' \item{fullalpha: }{Alpha estimated from the fastRUVIII model.}
#' \item{M: }{Replicate matrix.}
#' }
#' @export
#' @examples
#' ## Loading example data
#' data('example_sce', package = 'scMerge')
#' ## Previously computed stably expressed genes
#' data('segList_ensemblGeneID', package = 'scMerge')
#' ## Running an example data with minimal inputs
#' library(SingleCellExperiment)
#' exprsMat <- scMerge2(exprsMat = logcounts(example_sce),
#' batch = example_sce$batch,
#' ctl = segList_ensemblGeneID$mouse$mouse_scSEG)
#' assay(example_sce, "scMerge2") <- exprsMat$newY
#' 
#' 
#' example_sce = scater::runPCA(example_sce, exprs_values = 'scMerge2')                                       
#' scater::plotPCA(example_sce, colour_by = 'cellTypes', shape = 'batch')



scMerge2 <- function(exprsMat, 
                     batch, 
                     cellTypes = NULL, 
                     condition = NULL,
                     ctl = rownames(exprsMat), 
                     chosen.hvg = NULL,
                     ruvK = 20, 
                     use_bpparam = BiocParallel::SerialParam(),
                     use_bsparam = BiocSingular::RandomParam(),
                     use_bnparam = BiocNeighbors::AnnoyParam(),
                     pseudoBulk_fn = "create_pseudoBulk",
                     k_pseudoBulk = 30, 
                     k_celltype = 10,
                     exprsMat_counts = NULL,
                     cosineNorm = TRUE,
                     return_subset = FALSE,
                     return_subset_genes = NULL,
                     return_matrix = TRUE,
                     verbose = TRUE,
                     seed = 1) {
    
    
    set.seed(seed)
    
    
    # Input check
    
    ## In case there are complete zeroes in the rows or columns
    colsum_exprs <- DelayedMatrixStats::colSums2(exprsMat)
    rowsum_exprs <- DelayedMatrixStats::rowSums2(exprsMat)
    # if(any(rowsum_exprs == 0)){
    #     message("Automatically removing ",
    #             sum(rowsum_exprs == 0), " genes that are all zeroes in the data")
    #     exprsMat <- exprsMat[rowsum_exprs != 0, ]
    #     if (!is.null(exprsMat_counts)) {
    #         exprsMat_counts <- exprsMat_counts[rowsum_exprs != 0, ]
    #     }
    # }
    # 
    
    if(any(colsum_exprs == 0)){
        message("Automatically removing ", sum(colsum_exprs == 0), 
                " cells that are all zeroes in the data")
        exprsMat <- exprsMat[, colsum_exprs != 0]
        batch <- batch[colsum_exprs != 0]
        
        if (!is.null(cellTypes)) {
            cellTypes <- cellTypes[colsum_exprs != 0]
        }
        
        if (!is.null(condition)) {
            condition <- condition[colsum_exprs != 0]
        }
        
        if (!is.null(exprsMat_counts)) {
            exprsMat_counts <- exprsMat_counts[ , colsum_exprs != 0]
        }
        
    }
    
    
    
    .check_input_scMerge2(exprsMat, batch, 
                          cellTypes, condition, 
                          ctl, chosen.hvg, return_subset_genes,
                          exprsMat_counts)
    
    
    
    batch <- as.character(batch)
    if (!is.null(condition)) {
        condition <- as.character(condition)
    }
    batch_list <- sort(unique(batch))
    if (!is.null(condition)) {
        cond_mode <- TRUE
        condition_list <- sort(unique(condition))
        sample <- paste(batch, condition, sep = "_")
        sample_list <- sort(unique(sample))
        pseudobulk_sample <- sample
        pseudobulk_sample_list <- sample_list
        pseudobulk_batch <- batch
        pseudobulk_batch_list <- batch_list
    } else {
        cond_mode <- FALSE
        sample <- batch
        pseudobulk_sample <- batch
        pseudobulk_sample_list <- batch_list
        pseudobulk_batch <- batch
        pseudobulk_batch_list <- batch_list
    }
    
    
    
    
    meta_sample <- cbind(batch = batch, sample = sample, condition = condition)
    meta_sample <- unique(meta_sample)
    
    if (methods::is(exprsMat, "matrix")) {
        exprsMat <- methods::as(exprsMat, "CsparseMatrix")
    }
    
    
    if (is.null(chosen.hvg)) {
        chosen.hvg <- feature_selection(exprsMat, batch, use_bpparam = use_bpparam)
    }
    
    
    
    clust <- group_wtihin_batch(exprsMat = exprsMat,
                                cellTypes = cellTypes,
                                pseudobulk_batch_list = pseudobulk_batch_list,
                                pseudobulk_batch = pseudobulk_batch,
                                pseudobulk_sample_list = pseudobulk_sample_list,
                                pseudobulk_sample = pseudobulk_sample,
                                chosen.hvg = chosen.hvg,
                                k_celltype = k_celltype,
                                use_bpparam = use_bpparam,
                                use_bnparam = use_bnparam,
                                seed = seed,
                                verbose = verbose)
    
    
    if (cosineNorm) {
        if (verbose) {
            print("Normalising data")
        }
        exprsMat <- batchelor::cosineNorm(exprsMat, BPPARAM = use_bpparam)
    }
    
    
    if (verbose) {
        print("Constructing pseudo-bulk")
    }
    
    construct_bulk_res <- construct_pseudoBulk(exprsMat = exprsMat,
                                               exprsMat_counts = exprsMat_counts,
                                               pseudobulk_sample = pseudobulk_sample,
                                               pseudobulk_sample_list = pseudobulk_sample_list,
                                               clust = clust,
                                               cosineNorm = cosineNorm,
                                               pseudoBulk_fn = pseudoBulk_fn,
                                               k_pseudoBulk = k_pseudoBulk, 
                                               use_bpparam = use_bpparam,
                                               verbose = verbose) 
    
    
    
    
    bulkExprs <- construct_bulk_res$bulkExprs
    pseudo_batch <- construct_bulk_res$pseudo_batch
    bulk_clustering_res <- construct_bulk_res$bulk_clustering_res
    bulk_clustering_distProp <- construct_bulk_res$bulk_clustering_distProp
    
    
    
    if (verbose) {
        print("Identifying MNC using pseudo-bulk:")
    }
    
    if (cond_mode) {
        
        if (verbose) {
            print("condition_mode")
        }
        
        
        
        
        bulk_clustering_res_condition <- meta_sample[match(pseudobulk_sample_list, meta_sample[, "sample"]), "condition"]
        
        replicate_vector_condition <- lapply(condition_list, function(x) {
            idx <- bulk_clustering_res_condition == x
            bulk_clustering_subset <- bulk_clustering_res[idx]
            mnc_res <- findMNC(bulkExprs[chosen.hvg, names(unlist(bulk_clustering_subset))], 
                               bulk_clustering_subset, dist = "cor",
                               BPPARAM = use_bpparam, 
                               plot_igraph = verbose)
            
            replicate_vector <- mncReplicate(clustering_list = bulk_clustering_subset, 
                                             clustering_distProp = bulk_clustering_distProp[idx], 
                                             replicate_prop = 1, mnc_df = mnc_res)
            
        })
        
        replicate_vector_condition <- lapply(1:length(replicate_vector_condition), function(i) {
            res <- paste(condition_list[i], replicate_vector_condition[[i]], sep = "_")
            names(res) <- names(replicate_vector_condition[[i]])
            res
        })
        
        replicate_vector_condition <- unlist(replicate_vector_condition)
        replicate_vector <- replicate_vector_condition[colnames(bulkExprs)]
    }
    
    if (!cond_mode) {
        
        
        mnc_res <- findMNC(bulkExprs[chosen.hvg, ], 
                           bulk_clustering_res, dist = "cor",
                           BPPARAM = use_bpparam, 
                           plot_igraph = verbose)
        
        
        replicate_vector <- mncReplicate(clustering_list = bulk_clustering_res, 
                                         clustering_distProp = bulk_clustering_distProp, 
                                         replicate_prop = 1, mnc_df = mnc_res)
    }
    
    
    
    if (return_subset) {
        if (is.null(return_subset_genes)) {
            subset_genes <- chosen.hvg
        } else {
            subset_genes <- intersect(return_subset_genes, rownames(exprsMat)) 
        }
    } else {
        subset_genes <- NULL
    }
    
    
    if (verbose) {
        print("Running RUV")
    }
    
    
    
    res <- pseudoRUVIII(Y = t(exprsMat), 
                        Y_pseudo = t(bulkExprs), 
                        M = ruv::replicate.matrix(replicate_vector), 
                        ctl = rownames(exprsMat) %in% ctl,
                        k = ruvK, 
                        BSPARAM = use_bsparam,
                        return.info = TRUE,
                        subset = subset_genes,
                        normalised = return_matrix)
    
    
    if (return_matrix) {
        res$newY <- t(res$newY)
    }
    
    
    gc(reset = TRUE)
    
    return(res)
    
}




#' @title getAdjustedMat
#' @description Get Adjusted Matrix with scMerge2 parameter estimated 
#' @author Yingxin Lin
#' @param exprsMat A gene (row) by cell (column) matrix to be adjusted.
#' @param fullalpha A matrix indicates the estimated alpha returned by \code{scMerge2()}.
#' @param ctl A character vector of negative control. It should have a non-empty intersection with the rows of exprsMat.
#' @param adjusted_means A rowwise mean of the gene by cell matrix 
#' @param ruvK An integer indicates the number of unwanted variation factors that are removed, default is 20.
#' @param return_subset_genes An optional character vector of indicates the subset of genes will be adjusted.
#' 
#' @importFrom DelayedArray DelayedArray t
#' @importFrom DelayedMatrixStats colMeans2
#' 
#' @return Returns the adjusted matrix will be return.
#' 
#' @export
#' @examples
#' ## Loading example data
#' data('example_sce', package = 'scMerge')
#' ## Previously computed stably expressed genes
#' data('segList_ensemblGeneID', package = 'scMerge')
#' ## Running an example data with minimal inputs
#' library(SingleCellExperiment)
#' scMerge2_res <- scMerge2(exprsMat = logcounts(example_sce),
#' batch = example_sce$batch,
#' ctl = segList_ensemblGeneID$mouse$mouse_scSEG,
#' return_matrix = FALSE)
#' cosineNorm_mat <- batchelor::cosineNorm(logcounts(example_sce))
#' adjusted_means <- DelayedMatrixStats::rowMeans2(cosineNorm_mat)
#' newY <- getAdjustedMat(cosineNorm_mat, scMerge2_res$fullalpha,
#'               ctl = segList_ensemblGeneID$mouse$mouse_scSEG,
#'               ruvK = 20,
#'               adjusted_means = adjusted_means)
#' assay(example_sce, "scMerge2") <- newY
#' 
#' example_sce = scater::runPCA(example_sce, exprs_values = 'scMerge2')                                       
#' scater::plotPCA(example_sce, colour_by = 'cellTypes', shape = 'batch')





getAdjustedMat <- function(exprsMat, fullalpha, 
                           ctl = rownames(exprsMat), 
                           adjusted_means = NULL,
                           ruvK = 20,
                           return_subset_genes = NULL) {
    
    exprsMat <- DelayedArray(exprsMat)
    exprsMat <- DelayedArray::t(exprsMat)
    if (class(ctl) %in% c("character")) {
        ctl <- intersect(colnames(exprsMat), ctl)
        if (length(ctl) == 0) {
            stop("No provided ctl genes in the data")
        }
    }
    
    if (!is.null(return_subset_genes)) {
        return_subset_genes <- intersect(colnames(exprsMat), return_subset_genes)
        if (length(return_subset_genes) == 0) {
            stop("No provided return_subset_genes genes in the data")
        }
    }
    
    m <- nrow(exprsMat)

    if (is.null(adjusted_means)) {
        Y_stand <- sweep(exprsMat, 2, DelayedMatrixStats::colMeans2(exprsMat), "-")
    } else {
        Y_stand <- sweep(exprsMat, 2, adjusted_means, "-")
    }
    
    fullalpha <- DelayedArray(fullalpha)
    alpha <- fullalpha[seq_len(min(ruvK, nrow(fullalpha))), , drop = FALSE]
    W <- Y_stand[, ctl] %*% DelayedArray::t(alpha[, ctl]) %*% solve(alpha[, ctl] %*% DelayedArray::t(alpha[, ctl]))
    
    if (!is.null(return_subset_genes)) {
        newY <- exprsMat[, return_subset_genes] - W %*% alpha[, return_subset_genes]
    } else {
        newY <- exprsMat - W %*% alpha
    }
    
    newY <- t(newY)
    return(newY)
}

