#' scRUVIII
#'
#' perform scRUVIII
#'
#'
#' @author Yingxin Lin, Kevin Wang
#' @param Y The unnormalised SC data. A m by n matrix, where m is the number of observations and n is the number of features.
#' @param M The replicate mapping matrix.
#' The mapping matrix has m rows (one for each observation), and each column represents a set of replicates. The (i, j)-th entry of the mapping matrix is 1 if the i-th observation is in replicate set j, and 0 otherwise.
#' See ruv::RUVIII for more details.
#' @param ctl An index vector to specify the negative controls. Either a logical vector of length n or a vector of integers.
#' @param fullalpha Not used. Please ignore.
#' @param k The number of unwanted factors to remove. This is inherited from the ruvK argument from the scMerge::scMerge function.
#' @param return.info Additional information relating to the computation of normalised matrix. We recommend setting this to true.
#' @param cell_type An optional vector indicating the cell type information for each cell in the batch-combined matrix. If it is \code{NULL}, pseudo-replicate procedure will be run to identify cell type.
#' @param batch Batch information inherited from the scMerge::scMerge function.
#' @param return_all_RUV Whether to return extra information on the RUV function, inherited from the scMerge::scMerge function
#' @param fast_svd If \code{TRUE} we will use the randomized SVD algorithm to speed up computation. Otherwise, we will use the standard SVD built into R.
#' @param rsvd_prop If \code{fast_svd = TRUE}, then rsvd_prop is used to control for the accuracy of the randomized SVD computation.
#' If a lower value is used on a lower dimensional data (say < 1000 cell) will potentially yield a less accurate computed result but with a gain in speed.
#' The default of 0.1 tends to achieve a balance between speed and accuracy.
#' @export


# A function to perform loacation/scale adjustement to data as the input of RUVIII
# which also provides the option to select optimal RUVk according to the silhouette coefficient
scRUVIII <- function(Y = Y,
                     M = M, ctl = ctl,
                     fullalpha = NULL,
                     k = k, return.info = TRUE,
                     cell_type = NULL,
                     batch = NULL,
                     return_all_RUV = TRUE,
                     fast_svd = FALSE,
                     rsvd_prop = 0.1) {

  ## Transpose the data, since RUV assumes columns are genes.
  Y <- t(Y)


  # geneMeanMat<-matrix(rep(rowMeans(Y),ncol(Y)),ncol=ncol(Y))
  # geneSdMat<-matrix(rep(apply(Y-geneMeanMat,1,sd),ncol(Y)),ncol=ncol(Y))
  # normY<-(Y-geneMeanMat)/geneSdMat

  if (is.null(batch)) {
    warning("No batch info!")
    return(NULL)
  }

  ## Standardise the data
  scale_res <- standardize(Y, batch)
  normY <- scale_res$s.data
  geneSdMat <- sqrt(scale_res$stand.var) %*% t(rep(1, ncol(Y)))
  geneMeanMat <- scale_res$stand.mean
  ################
  ## We will always run an initial RUV, based on whether fast_svd is TRUE or not.

  ruv3_initial <- fastRUVIII(
    Y = t(normY),
    ctl = ctl,
    k = k[1],
    M = M,
    fullalpha = fullalpha,
    return.info = return.info,
    fast_svd = fast_svd,
    rsvd_prop = rsvd_prop
  )

  ruv3_initial$k <- k
  ## Finitial initial RUV3 run
  ###################
  ## The final computed result is always be ruv3res_list.
  ## If we have only one ruvK value, then the result is ruv3res_list, with only one element, corresponding to our initial run.
  ruv3res_list = vector("list", length = length(k))
  ruv3res_list[[1]] = ruv3_initial

  if(length(k) == 1){

  } else {
    ## If we have more than one ruvK value, then we feed the result to the ruv::RUVIII function
    ## (there is no need for fast_svd, since we already have the fullalpha)
    for(i in 2:length(k)){
      ruv3res_list[[i]] = fastRUVIII(
        Y = t(normY),
        ctl = ctl,
        k = k[i],
        M = M,
        fullalpha = ruv3_initial$fullalpha,
        return.info = return.info,
        fast_svd = FALSE
      )
    } ## End for loop
  } ## End else(length(k) == 1)

  names(ruv3res_list) = k
  ##################
  ## No need to run for length(k)==1
  if(length(k) == 1){
    f_score <- 1
    names(f_score) <- k
  } else {
    ## Cell type information

    cat("Selecting optimal RUVk\n")

    if (is.null(cell_type)) {
      cat("No cell type info, replicate matrix will be used as cell type info\n")
      cell_type <- apply(M, 1, function(x) which(x == 1))
    }
    ##################
    ## Computing the silhouette coefficient from kBET package
    sil_res <- do.call(cbind, lapply(ruv3res_list,
                                     FUN = function(x) {
                                       ## Computing the 10 PCA vectors using rsvd::rpca
                                       pca.data <- irlba::prcomp_irlba(x$newY, n = 10)
                                       # pca.data <- rsvd::rpca(x$newY, k = 10, rand = 1)
                                       # pca.data<-prcomp(x$newY)
                                       c(
                                         kBET::batch_sil(pca.data, as.numeric(as.factor(cell_type))),
                                         kBET::batch_sil(pca.data, as.numeric(as.factor(batch)), nPCs = 10)
                                       )
                                     }
    ))
    ##################
    ## Computing the F scores based on the 2 silhouette coefficients
    f_score <- rep(NA, ncol(sil_res))

    for (i in 1:length(k)) {
      f_score[i] <- f_measure(zeroOneScale(sil_res[1, ])[i], 1 - zeroOneScale(sil_res[2, ])[i])
    }
    names(f_score) <- k

    # print(sil_res)
    print(paste("optimal ruvK:", k[which.max(f_score)]))

    ## Not showing, if this needs displaying, consider implementing ggplot.
    plot(k, f_score, pch = 16, col = "light grey")
    lines(k, f_score)
    points(ruv3_initial$k, f_score[ruv3_initial$k], col = "red", pch = 16)

  }

  ##################
  for (i in 1:length(ruv3res_list)) {
    ruv3res_list[[i]]$newY <- t((t(ruv3res_list[[i]]$newY) * geneSdMat + geneMeanMat))
  }
  ##################
  ## ruv3res is the normalised matrix having the maximum F-score
  ruv3res_optimal <- ruv3res_list[[which.max(f_score)]]



  if (return_all_RUV) {
    ## If return_all_RUV is TRUE, we will un-scale every normalised matrices
    ruv3res_list$optimal_ruvK = k[which.max(f_score)] ## Always record the optimal k
    return(ruv3res_list)
  } else {
    ## If reurn_all is FALSE, we will un-scale the only normalised matrices
    ruv3res_optimal$optimal_ruvK <- k[which.max(f_score)] ## Always record the optimal k
    return(ruv3res_optimal)
  }
} ## End scRUVIII function








#######################################################
zeroOneScale <- function(v) {
  v <- (v+1)/2
  return(v)
}
#######################################################
standardize <- function(exprsMat, batch) {
  num_cell <- ncol(exprsMat)
  num_batch <- length(unique(batch))
  batch <- as.factor(batch)
  grand.mean <- matrix(rowMeans(exprsMat), nrow = 1)
  stand.mean <- t(grand.mean) %*% t(rep(1, num_cell))
  design <- model.matrix(~-1 + batch)
  B.hat <- solve(crossprod(design), tcrossprod(t(design), as.matrix(exprsMat)))
  # var.pooled <- ((exprsMat - t(design %*% B.hat))^2) %*% rep(1 / num_cell, num_cell)
  var.pooled <- ((exprsMat - t(design %*% B.hat))^2) %*% rep(1 / (num_cell-num_batch), num_cell)
  s.data <- (exprsMat - stand.mean) / (sqrt(var.pooled) %*% t(rep(1, num_cell)))
  return(res = list(s.data = s.data, stand.mean = stand.mean, stand.var = var.pooled))
}
#######################################################
f_measure <- function(celltypes, batch) {
  f <- 2 * (celltypes * batch) / (celltypes + batch)
  return(f)
}
