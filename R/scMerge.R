#' @title Perform the scMerge algorithm 
#' @description Merge single-cell RNA-seq data from different batches and experiments 
#' leveraging (pseudo)-replicates and control genes.
#'
#' @param sce_combine A \code{SingleCellExperiment} object contains the batch-combined matrix with batch info in colData.
#' @param ctl A character vector of negative control. It should have a non-empty intersection with the rows of sce_combine.
#' @param kmeansK A vector indicates the kmeans's K for each batch. The length of kmeansK needs to be the same as the number of batch.
#' @param exprs A string indicating the name of the assay requiring batch correction in sce_combine, default is logcounts.
#' @param hvg_exprs A string indicating the assay that to be used for highly variable genes identification in sce_combine, default is counts.
#' @param batch_name A character indicating the name of the batch column, default to "batch"
#' @param marker An optional vector of markers, to be used in calculation of mutual nearest cluster. If no markers input, highly variable genes will be used instead.
#' @param marker_list An optional list of markers for each batch, which will be used in calculation of mutual nearest cluster.
#' @param ruvK An optional integer/vector indicating the number of unwanted variation factors that are removed, default is 20.
#' @param replicate_prop A number indicating the ratio of cells that are included in pseudo-replicates, ranges from 0 to 1. Default to 1.
#' @param cell_type An optional vector indicating the cell type information for each cell in the batch-combined matrix. If it is \code{NULL}, pseudo-replicate procedure will be run to identify cell type.
#' @param cell_type_match An optional logical input for whether to find mutual nearest cluster using cell type information.
#' @param cell_type_inc An optional vector indicating the indices of the cells that will be used to supervise the pseudo-replicate procedure.
#' @param svd_k If BSPARAM is set to \code{RandomParam} or \code{IrlbaParam} class from \code{BiocSingular} package, then 
#' \code{svd_k} will be used to used to reduce the computational cost of singular value decomposition. Default to 50.
#' @param BSPARAM A \code{BiocSingularParam} class object from the \code{BiocSingular} package is used. Default is ExactParam().
#' @param dist The distance metrics that are used in the calculation of the mutual nearest cluster, default is Pearson correlation.
#' @param WV A optional vector indicating the wanted variation factor other than cell type info, such as cell stages.
#' @param WV_marker An optional vector indicating the markers of the wanted variation.
#' @param BPPARAM A \code{BiocParallelParam} class object from the \code{BiocParallel} package is used. Default is SerialParam().
#' @param return_all_RUV If \code{FALSE}, then only returns a \code{SingleCellExperiment} object with original data and one normalised matrix.
#' Otherwise, the \code{SingleCellExperiment} object will contain the original data and one normalised matrix for \code{each} ruvK value. In this latter case, assay_name must have the same length as ruvK.
#' @param assay_name The assay name(s) for the adjusted expression matrix(matrices). If \code{return_all_RUV = TRUE} assay_name must have the same length as ruvK.
#' @param plot_igraph If \code{TRUE}, then during the un/semi-supervised scMerge, igraph plot will be displayed
#' @param verbose If \code{TRUE}, then all intermediate steps will be shown. Default to \code{FALSE}.
#' @return Returns a \code{SingleCellExperiment} object with following components:
#' \itemize{
#' \item{assays: }{the original assays and also the normalised matrix}
#' \item{metadata: }{containing the ruvK vector, ruvK_optimal based on F-score, and the replicate matrix}
#' }
#' @author Yingxin Lin, Kevin Wang
#' @import SummarizedExperiment
#' @importFrom S4Vectors metadata
#' @importFrom BiocParallel SerialParam
#' @importFrom BiocParallel bpparam
#' @importFrom DelayedArray rowSums
#' @importFrom DelayedArray colSums
#' @importFrom BiocSingular ExactParam
#' @export
#' @examples
#' ## Loading example data
#' data('example_sce', package = 'scMerge')
#' ## Previously computed stably expressed genes
#' data('segList_ensemblGeneID', package = 'scMerge')
#' ## Running an example data with minimal inputs
#' sce_mESC <- scMerge(sce_combine = example_sce,
#' ctl = segList_ensemblGeneID$mouse$mouse_scSEG,
#' kmeansK = c(3, 3),
#' assay_name = 'scMerge')
#' 
#' sce_mESC = scater::runPCA(sce_mESC, exprs_values = "logcounts")                      
#' scater::plotPCA(sce_mESC, colour_by = 'cellTypes', shape = 'batch')
#' 
#' sce_mESC = scater::runPCA(sce_mESC, exprs_values = 'scMerge')                                       
#' scater::plotPCA(sce_mESC, colour_by = 'cellTypes', shape = 'batch')


scMerge <- function(sce_combine, ctl = NULL, kmeansK = NULL, 
    exprs = "logcounts", hvg_exprs = "counts", batch_name = "batch", marker = NULL, 
    marker_list = NULL, ruvK = 20, replicate_prop = 1, cell_type = NULL, 
    cell_type_match = FALSE, cell_type_inc = NULL, BSPARAM = ExactParam(), 
    svd_k = 50, dist = "cor", WV = NULL, WV_marker = NULL, 
    BPPARAM = SerialParam(), return_all_RUV = FALSE,
    assay_name = NULL, plot_igraph = TRUE, verbose = FALSE) {
    
    
    ## In case there are complete zeroes in the rows or columns
    colsum_exprs = DelayedMatrixStats::colSums2(SummarizedExperiment::assay(sce_combine, exprs))
    rowsum_exprs = DelayedMatrixStats::rowSums2(SummarizedExperiment::assay(sce_combine, exprs))
    if(any(colsum_exprs == 0) | any(rowsum_exprs == 0)){
        message("Automatically removing ", sum(colsum_exprs == 0), " cells and ",
                sum(rowsum_exprs == 0), " genes that are all zeroes in the data")
        sce_combine = sce_combine[rowsum_exprs != 0, colsum_exprs != 0]
    }
    
    ## This is a very niche option, consider deprecation in the future
    if (length(ruvK) > 1){
        message("You chose more than one ruvK. The argument return_all_RUV is forced to be TRUE.")
        return_all_RUV = TRUE
    }
    
    check_input(sce_combine = sce_combine, exprs = exprs, hvg_exprs = hvg_exprs, 
                assay_name = assay_name, batch_name = batch_name, ruvK = ruvK, 
                return_all_RUV = return_all_RUV, cell_type = cell_type)
    
    ## Extracting data matrix from SCE object
    exprs_mat <- SummarizedExperiment::assay(sce_combine, exprs)
    
    
    ## Checking negative controls input
    sce_rownames <- rownames(sce_combine)
    if (is.null(ctl)) {
        stop("Negative control genes are needed. \n 
             You could consider to use data(segList) as a pre-curated list.")
    } else {
        if (is.character(ctl)) {
            ctl <- which(sce_rownames %in% ctl)
        }
        if (length(ctl) == 0) {
            stop("Negative control genes are needed. \n 
             You could consider to use data(segList) as a pre-curated list.",
                call. = FALSE)
        }
    }
    
    ## Dealing with batch
    sce_combine$batch = sce_combine[[batch_name]]
    
    if (is.factor(sce_combine$batch)) {
        batch <- droplevels(sce_combine$batch)
    } else {
        batch <- sce_combine$batch
    }
    
    
    ## Finding pseudo-replicates
    t1 <- Sys.time()
    repMat <- scReplicate(sce_combine = sce_combine, batch = batch, 
        kmeansK = kmeansK, exprs = exprs, hvg_exprs = hvg_exprs, 
        marker = marker, marker_list = marker_list, replicate_prop = replicate_prop, 
        cell_type = cell_type, cell_type_match = cell_type_match, 
        cell_type_inc = cell_type_inc, dist = dist, WV = WV, 
        WV_marker = WV_marker, BPPARAM = BPPARAM, 
        BSPARAM = BSPARAM, plot_igraph = verbose, verbose = verbose)
    t2 <- Sys.time()
    
    timeReplicates <- t2 - t1
    
    cat("Dimension of the replicates mapping matrix: \n")
    print(dim(repMat))
    
    ## Performing RUV normalisation
    
    message("Step 2: Performing RUV normalisation. This will take minutes to hours. \n")
    
    ruv3res <- scRUVIII(Y = exprs_mat, M = repMat, ctl = ctl, 
        k = ruvK, batch = batch, fullalpha = NULL, cell_type = cell_type, 
        return_all_RUV = return_all_RUV, BSPARAM = BSPARAM, BPPARAM = BPPARAM,
        svd_k = svd_k)
    t3 <- Sys.time()
    
    timeRuv <- t3 - t2
    
    sce_combine <- sce_combine
    
    if (return_all_RUV) {
        ## if return_all_RUV is TRUE, then the previous check ensured
        ## assay_name is not NULL and matches the length of ruvK And
        ## the scRUVIII should've just returned with a single result
        ## (ruv3res_optimal)
        listNewY <- lapply(ruv3res[names(ruv3res) != "optimal_ruvK"], 
                           function(x) {t(x$newY)})
        
        for (i in seq_len(length(listNewY))) {
            SummarizedExperiment::assay(sce_combine, assay_name[i]) <- listNewY[[i]]
        }
    } else {
        ## If return_all_RUV is FALSE, then scRUVIII should've just
        ## returned with a single result (ruv3res_optimal)
        SummarizedExperiment::assay(sce_combine, assay_name) <- t(ruv3res$newY)
    }
    
    S4Vectors::metadata(sce_combine) <- c(S4Vectors::metadata(sce_combine), 
        list(ruvK = ruvK, ruvK_optimal = ruv3res$optimal_ruvK, 
            scRep_res = repMat, timeReplicates = timeReplicates, 
            timeRuv = timeRuv))
    
    message("scMerge complete!")
    
    return(sce_combine)
}  ## End scMerge function


check_input = function(sce_combine, exprs, hvg_exprs, assay_name, batch_name, ruvK, return_all_RUV, cell_type){
    #### Checking input expression
    if (is.null(exprs) | !exprs %in% SummarizedExperiment::assayNames(sce_combine)) {
        stop("The 'exprs' argument is NULL or it not a part of the supplied SingleCellExperiment object 'assayNames'")
    }
    
    if (is.null(hvg_exprs) | !hvg_exprs %in% SummarizedExperiment::assayNames(sce_combine)) {
        stop("The 'hvg_exprs' argument is NULL or it not a part of the supplied SingleCellExperiment object 'assayNames'")
    }
    
    #### Checking if the cell names are non-unique
    cell_names = colnames(sce_combine)
    
    if (length(cell_names) != length(unique(cell_names)) | is.null(cell_names)) {
        stop("Column names of the input SingleCellExperiment object must not contain duplicates nor NULL")
    }
    
    #### Returned assay name. Consider setting a default in future releases
    if (is.null(assay_name)) {
        stop("assay_name is NULL, please provide a name to store the results under")
    }
    
    
    #### Checking the batch info
    if (!(batch_name %in% colnames(colData(sce_combine)))) {
        stop(cat("Could not find a ", batch_name, " column in colData(sce_combine)"), 
             call. = FALSE)
    }
    
    if(sum(is.na(sce_combine[[batch_name]])) != 0){
        stop("NA's found the batch column, please remove")
    }
    
    
    #### Checking the cell_type info
    
    if(sum(is.na(cell_type)) != 0){
        stop("NA's found the cell_type input, please remove")
    }
    
    #### This is a very niche option, consider deprecation in the future
    if (return_all_RUV) {
        message("You chose return_all_RUV = TRUE. The result will contain all RUV computations. This could be a very large object.")
        ## We need an assay_name for every ruvK, if return_all_RUV is
        ## TRUE
        if (length(assay_name) != length(ruvK)) {
            stop("You chose return_all_RUV = TRUE. In this case, the length of assay_name must be equal to the length of ruvK")
        }
    }
}