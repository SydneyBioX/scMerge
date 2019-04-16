#' @title Simulate a simple matrix or SingleCellExperiment to test scMerge
#' @param m Number of observations
#' @param n Number of features
#' @param nc Number of negative controls
#' @param nCelltypes Number of cell-types
#' @param nBatch Number of batches
#' @param k Number of unwanted factors in simulation
#' @param lambda Rate parameter for random Poisson generation
#' @param sce If \code{TRUE}, returns a SingleCellExperiment object
#' @description This function is designed to generate Poisson-random-variable data matrix
#' to test on the internal algorithms of scMerge. It does not represent real biological situations 
#' and it is not intended to be used by end-users.
#' @return
#' If \code{sce} is FALSE, then the output is a list consists of
#' \itemize{
#' \item Y,  expression matrix generated through Poisson random variables,
#' \item ctl, a logical vector indicating the control genes,
#' \item M, replicate mapping matrix,
#' \item cellTypes, a vector indicating simulated cell types
#' \item batch, a vector indicating simulated batches
#' }
#' if \code{sce} is TRUE, a SingleCellExperiment wrapper will be applied on all above simulated objects.
#' @examples
#' set.seed(1)
#' L = ruvSimulate(m = 200, n = 1000, nc = 200, 
#' nCelltypes = 3, nBatch = 2, lambda = 0.1, k = 10, sce = TRUE)
#' print(L)
#' example <- scMerge(sce_combine = L,
#'                       ctl = paste0('gene', 1:500),
#'                       cell_type = L$cellTypes,
#'                       ruvK = 10,
#'                       assay_name = 'scMerge')
#'                       
#' scater::plotPCA(L, colour_by = 'cellTypes', shape = 'batch',
#'                  run_args = list(exprs_values = 'logcounts'))
#'                  
#' scater::plotPCA(example, colour_by = 'cellTypes', shape = 'batch',
#'                  run_args = list(exprs_values = 'scMerge'))
#' @export

ruvSimulate = function(m = 100, n = 5000, nc = floor(n/2), nCelltypes = 3, nBatch = 2, k = 20, lambda = 0.1, sce = FALSE) {
    ## m Number of observations n Number of features nc Number of negative controls
    p = 1
    ctl = rep(FALSE, n)
    ctl[seq_len(nc)] = TRUE
    celltypesNumbers = rep(seq_len(nCelltypes), length.out = m)
    X = matrix(celltypesNumbers, m, p)  ## Each number is an unique cell type
    # beta = matrix(stats::rpois(p * n, lambda = lambda), p, n)
    beta = matrix(stats::rpois(p * n, lambda = 2 * lambda), p, n)
    beta[, ctl] = 0
    
    
    W = matrix(stats::rpois(m * k, lambda = lambda), m, k)
    W[, 1] = matrix(seq_len(nBatch) - 1L, m, 1, byrow = FALSE)
    alpha = matrix(stats::rpois(k * n, lambda = lambda), k, n)
    
    
    epsilon = matrix(stats::rpois(m * n, lambda = lambda), m, n)
    
    Y = X %*% beta + W %*% alpha + epsilon
    
    rownames(Y) = paste0("cell", seq_len(m))
    colnames(Y) = paste0("gene", seq_len(n))
    
    cellTypes = paste0("cellTypes", rep_len(seq_len(nCelltypes), length.out = m))
    batch = paste0("batch", rep_len(seq_len(nBatch), length.out = m))
    
    M = ruv::replicate.matrix(data.frame(cellTypes, batch))
    
    if (sce) {
        result = SingleCellExperiment::SingleCellExperiment(assays = list(counts = t(Y), logcounts = log2(t(Y) + 1L)), colData = data.frame(cellTypes, 
            batch), metadata = list(M = M))
    } else {
        result = list(Y = Y, ctl = ctl, M = M, cellTypes = cellTypes, batch = batch)
    }
    return(result)
}

