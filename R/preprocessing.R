#' Combine SingleCellExperiment Object.
#'
#' Combind the \code{SingleCellExperiment} objects from different batches/experiments
#'
#' @param sce_list A list contains the \code{SingleCellExperiment} Object from each batch
#' @param method A string indicates the method of combining the gene expression matrix,
#' either \code{union} or \code{intersect}
#' @param cut_off_batch A numeric indicates cut-off for the proportion of a gene is expressed within each batch
#' @param cut_off_overall A numeric indicates the cut-off for the proportion of a gene is expressed overall
#' @param exprs A string vector indicates the exprssion matrix that are combined
#' @param colData_names A string vector indicates the \code{colData} that are combined
#' @param batch_names A string vector indicates the batch names
#'
#' @return A combined \code{SingleCellExperiment} objects
#' @author Yingxin Lin
#' @examples
#' require(SingleCellExperiment)
#' data("sce_mESC", package = "scMerge.data")
#' batchNames<-unique(sce_mESC$batch)
#' sce_list<-list(sce_mESC[,sce_mESC$batch==batchNames[1]],
#'                sce_mESC[,sce_mESC$batch==batchNames[2]],
#'                sce_mESC[,sce_mESC$batch==batchNames[3]],
#'                sce_mESC[,sce_mESC$batch==batchNames[4]],
#'                sce_mESC[,sce_mESC$batch==batchNames[5]])
#' sce_combine<-sce_cbind(sce_list,batch_names=batchNames)
#' @export


sce_cbind <- function(sce_list,
                      method = "union",
                      cut_off_batch = 0.01,
                      cut_off_overall = 0.01,
                      exprs = c("counts", "logcounts"),
                      colData_names = NULL,
                      batch_names = NULL) {
  n_batch <- length(sce_list)
  zero_list <- lapply(sce_list, function(x) apply(assay(x, exprs[1]), 1, function(z) sum(z == 0) / length(z)))
  expressed_list <- lapply(zero_list, function(x) names(which(x <= (1 - cut_off_batch))))
  for (i in 1:n_batch) {
    sce_list[[i]] <- sce_list[[i]][expressed_list[[i]], ]
  }

  if (is.null(method)) {
    method <- "union"
  }

  if (method == "intersect") {
    keep <- Reduce(intersect, expressed_list)
    sce_list <- lapply(sce_list, function(x) x[keep, ])
    assay_list <- list()
    for (i in 1:length(exprs)) {
      assay_list[[i]] <- do.call(cbind, lapply(sce_list, function(y) assay(y, exprs[i])))
    }
    names(assay_list) <- exprs
    colData_list <- do.call(rbind, lapply(sce_list, function(y) colData(y)[, colData_names, drop = F]))
    sce_combine <- SingleCellExperiment(
      assay = assay_list,
      colData = colData_list
    )
  } else {
    keep <- Reduce(union, expressed_list)

    assay_list <- list()
    for (i in 1:length(exprs)) {
      assay_list[[i]] <- do.call(cbind, lapply(sce_list, function(x) {
        mat <- matrix(0,
          nrow = length(keep),
          ncol = ncol(x),
          dimnames = list(
            keep,
            colnames(x)
          )
        )
        mat[rownames(x), ] <- assay(x, exprs[i])
        mat
      }))
    }
    names(assay_list) <- exprs
    colData_list <- do.call(rbind, lapply(sce_list, function(y) colData(y)[, colData_names, drop = F]))
    sce_combine <- SingleCellExperiment(
      assay = assay_list,
      colData = colData_list
    )

    zero_cbind <- apply(assay(sce_combine, exprs[1]), 1, function(z) sum(z == 0) / length(z))
    sce_combine <- sce_combine[names(zero_cbind[zero_cbind <= (1 - cut_off_overall)]), ]
  }

  if (is.null(batch_names)) {
    batch_names <- paste("batch", c(1:n_batch))
  }

  sce_combine$batch <- rep(batch_names, unlist(lapply(sce_list, ncol)))

  return(sce_combine)
}
