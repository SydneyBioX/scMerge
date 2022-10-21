#' @title Perform the scMerge2 (hierarchical) algorithm 
#' @description Merge single-cell data hierarchically from different batches and experiments 
#' leveraging (pseudo)-replicates, control genes and pseudo-bulk.
#' @author Yingxin Lin
#' @param exprsMat A gene (row) by cell (column) log-transformed matrix to be adjusted.
#' @param batch_list A list indicating the batch information for each cell in the batch-combined matrix.
#' @param h_idx_list A list indicating the indeces information in the hierarchical merging.
#' @param cellTypes An optional vector indicating the cell type information for each cell in the batch-combined matrix. 
#' If it is \code{NULL}, pseudo-replicate procedure will be run to identify cell type.
#' @param condition An optional vector indicating the condition information for each cell in the batch-combined matrix. 
#' @param ctl A character vector of negative control. It should have a non-empty intersection with the rows of exprsMat
#' @param chosen.hvg An optional character vector of highly variable genes chosen.
#' @param ruvK_list An integer indicates the number of unwanted variation factors that are removed, default is 20.
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
#' @export
#' @examples
#' 
#' ## Loading example data
#' data('example_sce', package = 'scMerge')
#' ## Previously computed stably expressed genes
#' data('segList_ensemblGeneID', package = 'scMerge')
#' 
#' 
#' # Create a fake sample information
#' example_sce$sample <- rep(c(1:4), each = 50)
#' 
#' # Construct a hierarchical index list
#' h_idx_list <- list(level1 = split(1:200, example_sce$batch),
#'                    level2 = list(1:200))
#' 
#' # Construct a batch information list
#' batch_list <- list(level1 = split(example_sce$sample, example_sce$batch),
#'                    level2 = list(example_sce$batch))
#' library(SingleCellExperiment)
#' exprsMat <- scMerge2h(exprsMat = logcounts(example_sce),
#' batch_list = batch_list,
#' h_idx_list = h_idx_list,
#' ctl = segList_ensemblGeneID$mouse$mouse_scSEG,
#' ruvK_list = c(2, 5))
#' assay(example_sce, "scMerge2") <- exprsMat[[length(h_idx_list)]]
#' example_sce = scater::runPCA(example_sce, exprs_values = 'scMerge2')                                       
#' scater::plotPCA(example_sce, colour_by = 'cellTypes', shape = 'batch')


scMerge2h <- function(exprsMat, 
                      batch_list = list(), 
                      h_idx_list = list(),
                      cellTypes = NULL, 
                      condition = NULL,
                      ctl = rownames(exprsMat), 
                      chosen.hvg = NULL,
                      ruvK_list = 20, 
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
    
    if (length(ruvK_list) == 1) {
        ruvK_list <- rep(ruvK_list, length(h_idx_list))
    }
    
    
    #
    
    
    # Some input check
    
    .check_input_scMerge2h(exprsMat, 
                           h_idx_list, 
                           batch_list, 
                           cellTypes, condition, 
                           ctl, chosen.hvg, return_subset_genes,
                           exprsMat_counts)
    
    # If only one level, run scMerge2() directly....
    
    if (length(h_idx_list) == 1) {
        
        warning("There is only one level in the input hierarchical index list, scMerge2() will be performed.")
        res <- scMerge2(exprsMat,
                        batch = batch, 
                        cellTypes = cellTypes, 
                        condition = condition,
                        ctl = ctl, 
                        chosen.hvg = chosen.hvg,
                        ruvK = ruvK_list[h_level], 
                        use_bpparam = use_bpparam,
                        use_bsparam = use_bsparam,
                        use_bnparam = use_bnparam,
                        pseudoBulk_fn = pseudoBulk_fn,
                        k_pseudoBulk = k_pseudoBulk, 
                        k_celltype = k_celltype,
                        exprsMat_counts = exprsMat_counts,
                        cosineNorm = cosineNorm,
                        return_subset = FALSE,
                        return_matrix = TRUE,
                        verbose = verbose,
                        seed = seed)
        return(res)
    }
    
    output_list <- list()
    for (h_level in seq_along(h_idx_list)) {
        
        output_list[[h_level]] <- list()
        
        
        for (h_data in seq_along(h_idx_list[[h_level]])) {
            
            if (verbose) {
                
                print(paste0("Hierarchical merging level ", h_level, ", data", h_data))
                
            }
            
            
            # if this is the first level integration
            if (h_level == 1) {
                
                idx <- h_idx_list[[h_level]][[h_data]]
                batch <- batch_list[[h_level]][[h_data]]
                
                #TODO: some input check for batch here
                res <- scMerge2(exprsMat[, idx],
                                batch = batch, 
                                cellTypes = cellTypes[idx], 
                                condition = condition[idx],
                                ctl = ctl, 
                                chosen.hvg = chosen.hvg,
                                ruvK = ruvK_list[h_level], 
                                use_bpparam = use_bpparam,
                                use_bsparam = use_bsparam,
                                use_bnparam = use_bnparam,
                                pseudoBulk_fn = pseudoBulk_fn,
                                k_pseudoBulk = k_pseudoBulk, 
                                k_celltype = k_celltype,
                                exprsMat_counts = exprsMat_counts[, idx],
                                cosineNorm = cosineNorm,
                                return_subset = FALSE,
                                return_matrix = TRUE,
                                verbose = verbose,
                                seed = seed)
                
                output_list[[h_level]][[h_data]] <- res$newY
                
                gc(reset = TRUE)
                
            } else {
                idx <- h_idx_list[[h_level]][[h_data]]
                batch <- batch_list[[h_level]][[h_data]]
                last_level_idx <- unlist(h_idx_list[[h_level - 1]])
                
                
                input_matrix <- list()
                last_level <- h_level - 1
                remain_idx <- idx
                while(last_level != 0) {
                    last_level_idx <- unlist(h_idx_list[[last_level]])
                    names(last_level_idx) <- NULL
                    common_idx <- intersect(last_level_idx, remain_idx)
                    
                    if (length(common_idx) != 0) {
                        input_matrix[[length(input_matrix) + 1]] <- output_list[[last_level]][, colnames(exprsMat)[common_idx]]
                    }
                    
                    remain_idx <- setdiff(idx, last_level_idx)
                    last_level <- last_level - 1
                }
                
                input_matrix <- do.call(DelayedArray::cbind, input_matrix)
                input_matrix <- cbind(input_matrix, 
                                      DelayedArray::DelayedArray(exprsMat[, colnames(exprsMat)[remain_idx]]))
                input_matrix <- input_matrix[, colnames(exprsMat)[idx]]
                
                # TODO: chosen.hvg selection
                output_list[[h_level]][[h_data]] <- scMerge2(input_matrix,
                                                             batch = batch, 
                                                             cellTypes = cellTypes[idx], 
                                                             condition = condition[idx],
                                                             ctl = ctl, 
                                                             chosen.hvg = rownames(input_matrix),
                                                             ruvK = ruvK_list[h_level], 
                                                             use_bpparam = use_bpparam,
                                                             use_bsparam = use_bsparam,
                                                             use_bnparam = use_bnparam,
                                                             pseudoBulk_fn = pseudoBulk_fn,
                                                             k_pseudoBulk = k_pseudoBulk, 
                                                             k_celltype = k_celltype,
                                                             exprsMat_counts = exprsMat_counts[, idx],
                                                             cosineNorm = cosineNorm,
                                                             return_subset = FALSE,
                                                             return_matrix = TRUE,
                                                             verbose = verbose,
                                                             seed = seed)$newY
                
                
            }
            
            
            
        }
        
        output_list[[h_level]] <- do.call(DelayedArray::cbind, output_list[[h_level]])
    }
    
    
    return(output_list)
    
}
