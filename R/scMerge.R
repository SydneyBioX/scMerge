#' Merge single-cell RNA-seq data from different batches, experiments, protocols.
#'
#' @param sce_combine A \code{SingleCellExperiment} object contains the batch-combined matrix with batch info in colData
#' @param ctl A vector indicates the indices of negative control
#' @param kmeansK A vector indicates the kmeans's K for each batch, length of KmeansK needs to be the same as the number of batch.
#' @param exprs A string inciates the assay that are used for batch correction, default is logcounts
#' @param hvg_exprs A string inciates the assay that are used for highly variable genes identification, default is counts
#' @param marker A vector of markers, which will be used in calcualtion of mutual nearest cluster. If no markers input, highly variable genes will be used instead
#' @param marker_list A list of markers for each batch, which will be used in calcualtion of mutual nearest cluster.
#' @param ruvK A vector indicates the number(s) of unwanted variation factors that are removed, default is 25.
#' @param replicate_prop A number indicates the ratio of cells that are included in pseudo-replicates, ranges from 0 to 1
#' @param cell_type A vector indicates the cell type information for each cell in the batch-combined matrix. If it is \code{NULL}, pseudo-replicate procedure will be run to identify cell type.
#' @param cell_type_match Whether find mutual nearest cluster using cell type information
#' @param cell_type_inc A vector indicates the indices of the cells that will be used to supervise the pseudo-replicate procedure
#' @param fast_svd If \code{TRUE}, rsvd will be used for singular value decomposition calculation.
#' @param dist The distance metrics that are used in the calcualtion of the mutual nearest cluster, default is Pearson correlation.
#' @param WV A vector indicates the wanted variation factor other than cell type info, such as cell stages.
#' @param WV_marker A vector indicates the markers of the wanted variation.
#' @param return_all_RUV whether return all the RUV results
#' @param assay_name The assay name for the adjusted expression matrix.
#' @param return_sce If \code{TRUE}, return the whole \code{SingleCellExperiment} object. If \code{FALSE}, return the adjusted expression matrix.
#'
#' @return If \code{return_sce} is \code{TRUE}, return the \code{SingleCellExperiment} object with the adjusted assay.
#'
#' If \code{return_sce} is \code{FALSE}, return the followings
#' \item{normalised_matrix }{adjusted expression matrix}
#' \item{scRUV }{fastRUVIII output}
#' \item{mnc }{mutual nearest clusters}
#' \item{kmeans }{kmeans results}
#' @author Yingxin Lin
#' @examples
#' require(SingleCellExperiment)
#' #Loading example data
#' data("sce_mESC")
#' data("segList_ensemblGeneID")
#' #scMerge
#' sce_mESC<- scMerge(sce_mESC,
#'                    ctl = segList_ensemblGeneID$mouse$mouse_scSEG,
#'                    kmeansK = c(1,3,3,1,1),
#'                    assay_name = "scRUV")
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
                    fast_svd = TRUE,
                    dist = "cor",
                    WV = NULL,
                    WV_marker = NULL,
                    return_all_RUV = FALSE,
                    assay_name = NULL,
                    return_sce = TRUE) {

  ## Checking input expression
  if(is.null(exprs)){
    stop("exprs is NULL.")
  }

  ## Checking input expression assay name in SCE object
  if(!exprs %in% assayNames(sce_combine)){
    stop(paste("No assay named", exprs))
  }

  ## Extracting data matrix from SCE object
  exprs_mat <- assay(sce_combine, exprs)
  sce_rownames <- rownames(sce_combine)

  ## Checking if any rows or columns are purely zeroes
  if (sum(rowSums(exprs_mat) == 0) != 0 | sum(colSums(exprs_mat) == 0) != 0) {
    stop("There are rows or columns that are all zeros. Please remove it.")
  }

  ## Checking negative controls input
  if (is.null(ctl)) {
    warning("No negative control genes \n")
    if (return_sce) {
      return(sce_combine)
    } else {
      return(NULL)
    }
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
  }else{
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


  ## Performing RUV normalisation

  cat("Performing RUV normalisation... ;) \n")

  ruv3res <- scRUVIII(Y = t(exprs_mat),
                              M = repMat,
                              ctl = ctl,
                              k = ruvK,
                              batch = sce_batch,
                              fullalpha = NULL,
                              cell_type = cell_type,
                              return.info = T,
                              return_all = return_all_RUV,
                              fast_svd = fast_svd)


  # YXL: NEED FIX
  ## Return results
  if (return_all_RUV) {
    return(res = list(
      scRUV = ruv3res,
      mnc = mnc,
      kmeans = kmeans_res
    ))
  }

  if (return_sce) {
    if (is.null(assay_name)) {
      assay_name <- paste("scMerge_RUVk", ruv3res$k, sep = "")
    }
    assay(sce_combine, assay_name) <- t(ruv3res$newY)
    cat(paste("Return assay named \n", assay_name))
    return(sce_combine)
  } else {
    return(res = list(
      normalised_matrix = ruv3res$newY,
      scRUV = ruv3res,
      mnc = mnc,
      kmeans = kmeans_res
    ))
  }
}
