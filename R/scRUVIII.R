#' @title scRUVIII: RUVIII algorithm optimised for single cell data
#'
#' @description A function to perform location/scale adjustment to data as the input of
#' RUVIII which also provides the option to select optimal RUVk according to the
#' silhouette coefficient
#'
#'
#' @author Yingxin Lin, Kevin Wang
#' @param Y The unnormalised SC data. A m by n matrix, where m is the number of observations and n is the number of features.
#' @param M The replicate mapping matrix.
#' The mapping matrix has m rows (one for each observation), and each column represents a set of replicates.
#' The (i, j)-th entry of the mapping matrix is 1 if the i-th observation is in replicate set j, and 0 otherwise.
#' See ruv::RUVIII for more details.
#' @param ctl An index vector to specify the negative controls.
#' Either a logical vector of length n or a vector of integers.
#' @param fullalpha Not used. Please ignore.
#' @param k The number of unwanted factors to remove. This is inherited from the ruvK argument from the scMerge::scMerge function.
#' @param cell_type An optional vector indicating the cell type information for each cell
#' in the batch-combined matrix. If it is \code{NULL},
#' pseudo-replicate procedure will be run to identify cell type.
#' @param batch Batch information inherited from the scMerge::scMerge function.
#' @param return_all_RUV Whether to return extra information on the RUV function, inherited from the scMerge::scMerge function
#' @param fast_svd If \code{TRUE}, fast algorithms will be used for singular value decomposition calculation via the \code{irlba} and \code{rsvd} packages.
#' We recommend using this option when the number of cells is large (e.g. more than 1000 cells).
#' @param rsvd_prop If \code{fast_svd = TRUE}, then \code{rsvd_prop} will be used to used to
#' reduce the computational cost of randomised singular value decomposition.
#' We recommend setting this number to less than 0.25 to achieve a balance between numerical accuracy and computational costs.
#' If a lower value is used on a lower dimensional data (say < 1000 cell) will potentially yield a
#' less accurate computed result but with a gain in speed.
#' The default of 0.1 tends to achieve a balance between speed and accuracy.
#' @return A list consists of:
#' \itemize{
#' \item{RUV-normalised matrices:} If k has multiple values, then the RUV-normalised matrices using
#' all the supplied k values will be returned.
#' \item{optimal_ruvK:} The optimal RUV k value as determined by silhouette coefficient.
#' }
#' @export
#' @examples
#' L = ruvSimulate(m = 200, n = 1000, nc = 100, nCelltypes = 3, nBatch = 2, lambda = 0.1, sce = FALSE)
#' Y = t(log2(L$Y + 1L)); M = L$M; ctl = L$ctl; batch = L$batch;
#' res = scRUVIII(Y = Y, M = M, ctl = ctl, k = c(5, 10, 15, 20), batch = batch)

scRUVIII <- function(Y = Y, M = M, ctl = ctl, fullalpha = NULL, 
                     k = k, cell_type = NULL, batch = NULL, return_all_RUV = TRUE, 
                     fast_svd = FALSE, rsvd_prop = 0.1) {
    
    if (is.null(batch)) {
        warning("No batch info!")
        return(NULL)
    }
    
    ## Standardise the data
    scale_res <- standardize2(Y, batch)
    normY <- scale_res$s.data
    geneSdMat <- sqrt(scale_res$stand.var) %*% t(rep(1, ncol(Y)))
                               
    # geneSdVec <- sqrt(scale_res$stand.var)
    geneMeanVec <- scale_res$stand.mean
    # geneMeanMat <- scale_res$stand.mean  %*% t(rep(1, ncol(Y)))
    ## We will always run an initial RUV, based on whether
    ## fast_svd is TRUE or not.
    
    ruv3_initial <- fastRUVIII(Y = t(normY), ctl = ctl, k = k[1], 
                               M = M, fullalpha = fullalpha, return.info = TRUE, fast_svd = fast_svd, 
                               rsvd_prop = rsvd_prop)
    
    ruv3_initial$k <- k
    ## The computed result is ruv3res_list.  If we have only one
    ## ruvK value, then the result is ruv3res_list with only one
    ## element, corresponding to our initial run.
    ruv3res_list = vector("list", length = length(k))
    ruv3res_list[[1]] = ruv3_initial
    
    if (length(k) == 1) {
        
    } else {
        ## If we have more than one ruvK value then we feed the result
        ## to the ruv::RUVIII function (there is no need for fast_svd,
        ## since we already have the fullalpha)
        for (i in 2:length(k)) {
            ruv3res_list[[i]] = fastRUVIII(Y = t(normY), ctl = ctl, 
                                           k = k[i], M = M, fullalpha = ruv3_initial$fullalpha, 
                                           return.info = TRUE, fast_svd = FALSE)
        }  ## End for loop
    }  ## End else(length(k) == 1)
    
    names(ruv3res_list) = k
    ## Caculate sil. coef and F-score to select the best RUVk
    ## value.  No need to run for length(k)==1
    if (length(k) == 1) {
        f_score <- 1
        names(f_score) <- k
    } else {
        ## Cell type information will be used for calculating sil.coef
        ## and F-score
        cat("Selecting optimal RUVk \n")
        
        if (is.null(cell_type)) {
            cat("No cell type info, replicate matrix will be used as cell type info \n")
            cell_type <- apply(M, 1, function(x) which(x == 1))
        }
        ## Computing the silhouette coefficient from kBET package
        sil_res <- do.call(cbind, lapply(ruv3res_list, FUN = calculateSil, 
                                         fast_svd = fast_svd, cell_type = cell_type, batch = batch))
        ## Computing the F scores based on the 2 silhouette
        ## coefficients
        f_score <- rep(NA, ncol(sil_res))
        
        for (i in seq_len(length(k))) {
            f_score[i] <- f_measure(zeroOneScale(sil_res[1, ])[i], 
                                    1 - zeroOneScale(sil_res[2, ])[i])
        }
        names(f_score) <- k
        
        message("optimal ruvK:", k[which.max(f_score)])
        
        ## Not showing
        graphics::plot(k, f_score, pch = 16, col = "light grey")
        graphics::lines(k, f_score)
        graphics::points(ruv3_initial$k[which.max(f_score)], 
                         f_score[[which.max(f_score)]], col = "red", pch = 16)
        
    }
    
    ## Add back the mean and sd to the normalised data
    for (i in seq_len(length(ruv3res_list))) {
        ruv3res_list[[i]]$newY <- t((t(ruv3res_list[[i]]$newY) *
                                         geneSdMat + geneMeanVec))
    }
    ## ruv3res_list is all the normalised matrices ruv3res_optimal
    ## is the one matrix having the maximum F-score
    ruv3res_optimal <- ruv3res_list[[which.max(f_score)]]
    
    
    
    if (return_all_RUV) {
        ## If return_all_RUV is TRUE, we will return all the
        ## normalised matrices
        ruv3res_list$optimal_ruvK <- k[which.max(f_score)]
        return(ruv3res_list)
    } else {
        ## If return_all_RUV is FALSE, we will return the F-score
        ## optimal matrix
        ruv3res_optimal$optimal_ruvK <- k[which.max(f_score)]
        return(ruv3res_optimal)
    }
}  ## End scRUVIII function








####################################################### 
zeroOneScale <- function(v) {
    v <- (v + 1)/2
    return(v)
}
####################################################### 
standardize <- function(exprsMat, batch) {
    exprsMat <- as.matrix(exprsMat)
    num_cell <- ncol(exprsMat)
    num_batch <- length(unique(batch))
    batch <- as.factor(batch)
    grand.mean <- matrix(base::rowMeans(exprsMat), nrow = 1)
    stand.mean <- t(grand.mean) %*% t(rep(1, num_cell))
    design <- stats::model.matrix(~-1 + batch)
    B.hat <- solve(crossprod(design), tcrossprod(t(design), as.matrix(exprsMat)))
    var.pooled <- ((exprsMat - t(design %*% B.hat))^2) %*% rep(1/(num_cell - 
                                                                      num_batch), num_cell)
    s.data <- (exprsMat - stand.mean)/(sqrt(var.pooled) %*% t(rep(1, 
                                                                  num_cell)))
    return(res = list(s.data = s.data, stand.mean = stand.mean, 
                      stand.var = var.pooled))
}
###############################
standardize2 <- function(Y, batch) {
    num_cell <- ncol(Y)
    num_batch <- length(unique(batch))
    batch <- as.factor(batch)
    stand.mean <- DelayedMatrixStats::rowMeans2(Y)
    design <- stats::model.matrix(~-1 + batch)
    B.hat <- solve(t(design) %*% design, 
                   t(Y %*% design))
    
    var.pooled <- matrix(DelayedMatrixStats::rowSums2(
        ((Y - t(B.hat) %*% t(design))^2)
    )/(num_cell - num_batch),
    ncol = 1)
    s.data.dem = sqrt(var.pooled) %*% matrix(1, nrow = 1, ncol = num_cell)
    s.data <- (Y - stand.mean)/s.data.dem
    return(res = list(s.data = s.data, stand.mean = stand.mean, 
                      stand.var = var.pooled))
}
####################################################### 
f_measure <- function(cell_type, batch) {
    f <- 2 * (cell_type * batch)/(cell_type + batch)
    return(f)
}
####################################################### 
calculateSil <- function(x, fast_svd, cell_type, batch) {
    if (fast_svd & !any(dim(x$newY) < 50)) {
        pca.data <- irlba::prcomp_irlba(x$newY, n = 10)
    } else {
        pca.data <- stats::prcomp(x$newY)
    }
    
    result = c(kBET_batch_sil(pca.data, as.numeric(as.factor(cell_type)), 
                              nPCs = 10), kBET_batch_sil(pca.data, as.numeric(as.factor(batch)), 
                                                         nPCs = 10))
    return(result)
}

kBET_batch_sil <- function(pca.data, batch, nPCs = 10) {
    ## This function was copied from kBET, which cannot be
    ## imported because it is only a GitHub package
    dd <- as.matrix(stats::dist(pca.data$x[, seq_len(nPCs)]))
    score_sil <- summary(cluster::silhouette(as.numeric(batch), 
                                             dd))$avg.width
    return(score_sil)
}

