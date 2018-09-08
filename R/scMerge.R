#' Merge single-cell RNA-seq data from different batches, experiments, protocols.
#'
#' @param sce_combine A \code{SingleCellExperiment} object contains the batch-combined matrix with batch info in colData.
#' @param ctl A chatacter vector of negative control. It should have a non-empty intersection with the rows of sce_combine.
#' @param kmeansK A vector indicates the kmeans's K for each batch. The length of kmeansK needs to be the same as the number of batch.
#' @param exprs A string inciating the name of the assay requiring batch correction in sce_combine, default is logcounts.
#' @param hvg_exprs A string inciating the assay that to be used for highly variable genes identification in sce_combine, default is counts.
#' @param marker An optional vector of markers, to be used in calculation of mutual nearest cluster. If no markers input, highly variable genes will be used instead.
#' @param marker_list An optional list of markers for each batch, which will be used in calculation of mutual nearest cluster.
#' @param ruvK An optional integer/vector indicating the number of unwanted variation factors that are removed, default is 20.
#' @param replicate_prop A number indicating the ratio of cells that are included in pseudo-replicates, ranges from 0 to 1.
#' @param cell_type An optional vector indicating the cell type information for each cell in the batch-combined matrix. If it is \code{NULL}, pseudo-replicate procedure will be run to identify cell type.
#' @param cell_type_match An optional logical input for whether to find mutual nearest cluster using cell type information.
#' @param cell_type_inc An optional vector indicating the indices of the cells that will be used to supervise the pseudo-replicate procedure.
#' @param fast_svd If \code{TRUE}, randomised singular value decomposition will be used for singular value decomposition calculation. We recommend using this option when the number of cells is large (e.g. > 1000).
#' @param rsvd_prop If \code{fast_svd = TRUE}, then rsvd_prop will be used to used to reduce the computational cost of randomised singular value decomposition. We recommend setting this number to less than 0.25 to achieve a balance between numerical accuracy and computational costs.
#' @param dist The distance metrics that are used in the calculation of the mutual nearest cluster, default is Pearson correlation.
#' @param WV A optional vector indicating the wanted variation factor other than cell type info, such as cell stages.
#' @param WV_marker An optional vector indicating the markers of the wanted variation.
#' @param return_all_RUV If \code{FALSE}, then only returns a \code{SingleCellExperiment} object with original data and one normalised matrix.
#' Otherwise, the \code{SingleCellExperiment} object will contain the original data and one normalised matrix for \code{each} ruvK value. In this latter case, assay_name must have the same length as ruvK.
#' @param assay_name The assay name(s) for the adjusted expression matrix(matrices). If \code{return_all_RUV = TRUE} assay_name must have the same length as ruvK.
#'
#' @return Returns a \code{SingleCellExperiment} object with following:
#'
#' \item{metadata}{containing the ruvK vector, ruvK_optimal based on F-score, and the replicate matrix}
#' \item{assays}{the original matrices and also the normalised matrices}
#' @author Yingxin Lin, Kevin Wang
#' @examples
#' suppressPackageStartupMessages({
#' library(SingleCellExperiment)
#' library(scater)
#' library(scMerge)
#' library(scMerge.data)
#' })
#' # Loading example data
#' data("sce_mESC", package = "scMerge.data")
#' # Previously computed stably expressed genes
#' data("segList_ensemblGeneID")
#' # Running an example data with minimal inputs
#' sce_mESC <- scMerge(
#'                       sce_combine = sce_mESC,
#'                       ctl = segList_ensemblGeneID$mouse$mouse_scSEG,
#'                       kmeansK = c(1,3,3,1,1),
#'                       assay_name = "scMerge")
#' scater::plotPCA(sce_mESC, colour_by = "cellTypes", shape = "batch",
#'                  run_args = list(exprs_values = "logcounts"))
#' scater::plotPCA(sce_mESC, colour_by = "cellTypes", shape = "batch",
#'                  run_args = list(exprs_values = "scMerge"))
#' @export



scMerge <- function(sce_combine,
                    ctl = NULL,
                    kmeansK = NULL,
                    exprs = "logcounts",
                    hvg_exprs = "counts",
                    marker = NULL,
                    marker_list = NULL,
                    ruvK = 20,
                    replicate_prop = 0.5,
                    cell_type = NULL,
                    cell_type_match = FALSE,
                    cell_type_inc = NULL,
                    fast_svd = FALSE,
                    rsvd_prop = 0.1,
                    dist = "cor",
                    WV = NULL,
                    WV_marker = NULL,
                    return_all_RUV = FALSE,
                    assay_name = NULL) {

  ## Checking input expression
  if(is.null(exprs)){
    stop("exprs is NULL.")
  }

  if(is.null(assay_name)){
    stop("assay_name is NULL")
  }

  if(return_all_RUV){
    message("You chose return_all_RUV = TRUE, the result will contain all RUV computations. This could be a very large object.")
    ## We need an assay_name for every ruvK, if return_all_RUV is TRUE
    if(length(assay_name) != length(ruvK)){
      stop("You chose return_all_RUV = TRUE. In this case, the length of assay_name must be equal to the length of ruvK")
    }
  }




  ## Checking input expression assay name in SCE object
  if(!exprs %in% assayNames(sce_combine)){
    stop(paste("No assay named", exprs))
  }

  ## Extracting data matrix from SCE object
  exprs_mat <- SummarizedExperiment::assay(sce_combine, exprs)
  sce_rownames <- rownames(sce_combine)

  ## Checking if any rows or columns are purely zeroes
  if (sum(rowSums(exprs_mat) == 0) != 0 | sum(colSums(exprs_mat) == 0) != 0) {
    stop("There are rows or columns that are all zeros. Please remove it.")
  }

  ## Checking negative controls input
  if (is.null(ctl)) {
    stop("Negative control genes are needed. \n")
  } else {
    if (class(ctl) == "character") {
      ctl <- which(sce_rownames %in% ctl)
    }
    if (length(ctl) == 0) {
      stop("Could not find any negative control genes in the row names of the expression matrix", call. = FALSE)
    }
  }

  ## Checking the batch info
  if (is.null(sce_combine$batch)) {
    stop("Could not find a 'batch' column in colData(sce_combine)", call. = FALSE)
  }

  if(class(sce_combine$batch) == "factor"){
    sce_batch <- droplevels(sce_combine$batch)
  } else {
    sce_batch <- sce_combine$batch
  }

  ## Finding pseudo-replicates
  repMat <- scReplicate(sce = sce_combine,
                        batch = sce_batch,
                        kmeansK = kmeansK,
                        exprs = exprs,
                        hvg_exprs = hvg_exprs,
                        marker = marker,
                        marker_list = marker_list,
                        replicate_prop = replicate_prop,
                        cell_type = cell_type,
                        cell_type_match = cell_type_match,
                        cell_type_inc = cell_type_inc,
                        dist = dist,
                        WV = WV,
                        WV_marker = WV_marker)


  cat("Dimension of the replicates mapping matrix \n")
  print(dim(repMat))
  cat("\n")

  ## Performing RUV normalisation

  cat("Performing RUV normalisation... This might take a few minutes... \n")

  ruv3res <- scRUVIII(Y = t(exprs_mat),
                      M = repMat,
                      ctl = ctl,
                      k = ruvK,
                      batch = sce_batch,
                      fullalpha = NULL,
                      cell_type = cell_type,
                      return.info = TRUE,
                      return_all_RUV = return_all_RUV,
                      fast_svd = fast_svd,
                      rsvd_prop = rsvd_prop)


  sce_final_result = sce_combine
  if(!return_all_RUV){
    ## If return_all_RUV is FALSE, then scRUVIII should've just returned with a single result (ruv3res_optimal)
    assay(sce_final_result, assay_name) <- t(ruv3res$newY)
  } else{
    ## if return_all_RUV is TRUE, then the previous check ensured assay_name is not NULL and matches the length of ruvK
    ## And the scRUVIII should've just returned with a single result (ruv3res_optimal)
    listNewY = lapply(ruv3res[names(ruv3res) != "optimal_ruvK"], function(x){t(x$newY)})
    for(i in 1:length(listNewY)){
      assay(sce_final_result, assay_name[i]) <- listNewY[[i]]
    } ## End for loop
  }

  metadata(sce_final_result) = c(
    metadata(sce_combine),
    list(
      "ruvK" = ruvK,
      "ruvK_optimal" = ruv3res$optimal_ruvK,
      "scRep_res" = repMat
    ))

  return(sce_final_result)
} ## End scMerge function
