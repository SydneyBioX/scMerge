#' @title Combind several \code{SingleCellExperiment} objects from different batches/experiments
#'
#' @description Combind several \code{SingleCellExperiment} objects from different batches/experiments.
#'
#' @param sce_list A list contains the \code{SingleCellExperiment} Object from each batch
#' @param method A string indicates the method of combining the gene expression matrix,
#' either \code{union} or \code{intersect}
#' @param cut_off_batch A numeric vector indicating the cut-off for the proportion of a gene is expressed within each batch
#' @param cut_off_overall A numeric vector  indicating the cut-off for the proportion of a gene is expressed overall data
#' @param exprs A string vector indicating the expression matrices to be combined. The first assay named will be used to determine the proportion of zeores.
#' @param colData_names A string vector indicating the \code{colData} that are combined
#' @param batch_names A string vector indicating the batch names for the output sce object
#' @import SingleCellExperiment
#'
#' @return A \code{SingleCellExperiment} object with the list of SCE objects combined.
#' @author Yingxin Lin
#' @export
#' @examples
#' library(SingleCellExperiment)
#' data('example_sce', package = 'scMerge')
#' batch_names<-unique(example_sce$batch)
#' sce_list<-list(example_sce[,example_sce$batch=='batch2'],
#'                example_sce[,example_sce$batch=='batch3'])
#' sce_combine<-sce_cbind(sce_list,batch_names=batch_names)


sce_cbind <- function(sce_list, method = NULL, cut_off_batch = 0.01, cut_off_overall = 0.01, 
    exprs = c("counts", "logcounts"), colData_names = NULL, batch_names = NULL) {
    
    message("The assay named '", exprs[1], "' will be used to determine the proportion of zeroes for each batch")
    
    n_batch <- length(sce_list)
    zero_list <- lapply(sce_list, function(x) apply(assay(x, exprs[1]), 1, function(z) mean(z == 
        0)))
    
    expressed_list <- lapply(zero_list, function(x) names(which(x <= (1 - cut_off_batch))))
    for (i in seq_len(n_batch)) {
        sce_list[[i]] <- sce_list[[i]][expressed_list[[i]], ]
    }
    
    
    
    if (is.null(method)) {
        message("As method is NULL, we will use the union method. \n")
        method = "union"
    }
    
    if (!method %in% c("union", "intersect")) {
        stop("The method parameter can only be either 'union' or 'intersect'. \n ")
    }
    
    if (method == "intersect") {
        keep <- Reduce(intersect, expressed_list)
        sce_list <- lapply(sce_list, function(x) x[keep, ])
        assay_list <- list()
        for (i in seq_len(length(exprs))) {
            assay_list[[i]] <- do.call(cbind, lapply(sce_list, function(y) assay(y, 
                exprs[i])))
        }
        names(assay_list) <- exprs
        colData_list <- do.call(rbind, lapply(sce_list, function(y) colData(y)[, colData_names, 
            drop = FALSE]))
        sce_combine <- SingleCellExperiment(assay = assay_list, colData = colData_list)
    }
    
    if (method == "union") {
        keep <- Reduce(union, expressed_list)
        
        assay_list <- list()
        for (i in seq_len(length(exprs))) {
            assay_list[[i]] <- do.call(cbind, lapply(sce_list, function(x) {
                mat <- matrix(0, nrow = length(keep), ncol = ncol(x), dimnames = list(keep, 
                  colnames(x)))
                mat[rownames(x), ] <- assay(x, exprs[i])
                return(mat)
            }))
        }
        names(assay_list) <- exprs
        colData_list <- do.call(rbind, lapply(sce_list, function(y) colData(y)[, colData_names, 
            drop = FALSE]))
        
        sce_combine <- SingleCellExperiment(assay = assay_list, colData = colData_list)
        
        zero_cbind <- apply(assay(sce_combine, exprs[1]), 1, function(z) mean(z == 
            0))
        sce_combine <- sce_combine[names(zero_cbind[zero_cbind <= (1 - cut_off_overall)]), 
            ]
    }
    
    if (is.null(batch_names)) {
        batch_names <- paste("batch", seq_len(n_batch))
    }
    
    sce_combine$batch <- rep(batch_names, unlist(lapply(sce_list, ncol)))
    
    return(sce_combine)
}
