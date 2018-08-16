#' scRUVIII
#'
#' perform scRUVIII
#'
#'
#' @author Yingxin Lin
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
#' @param return_all Whether to return extra information on the RUV function, inherited from the scMerge::scMerge function
#' @param fast_svd If \code{TRUE} we will use the randomized SVD algorithm to speed up computation. Otherwise, we will use the standard SVD built into R.
#' @param propEigen If \code{fast_svd = TRUE}, then propEigen is used to control for the accuracy of the randomized SVD computation.
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
                     return_all = F,
                     fast_svd = FALSE,
                     propEigen = 0.1) {

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



  if (length(k) == 1) { ## If user only specified one RUVk value
    if (fast_svd) { ## And the user wanted the rsvd option, we will provide that
      ruv3res <- fastRUVIII(
        Y = t(normY),
        ctl = ctl,
        k = k,
        M = M,
        fullalpha = fullalpha,
        return.info = return.info,
        propEigen = propEigen
      )
    } else { ## If the user wanted the usual svd option, we will provide that here
      ruv3res <- ruv::RUVIII(
        Y = t(normY),
        ctl = ctl,
        k = k,
        M = M,
        fullalpha = fullalpha,
        return.info = return.info
      )
    }

    ruv3res$k <- k
  } else { ## If user only specified more than one RUVk value, we will perform a computation with a single value first
    ruv3res_list <- list()
    print(paste("k =", k[1]))
    if (fast_svd) {
      ruv3res_list[[1]] <- fastRUVIII(
        Y = t(normY),
        ctl = ctl,
        k = k[1],
        M = M,
        fullalpha = fullalpha,
        return.info = return.info,
        propEigen = propEigen
      )
    } else {
      ruv3res_list[[1]] <- RUVIII(
        Y = t(normY),
        ctl = ctl,
        k = k[1],
        M = M,
        fullalpha = fullalpha,
        return.info = return.info
      )
    } ## End fast_svd criterion
    ## If user only specified more than one RUVk value, then, all subsequent normalisation can be calculated using the first decomposition value.
    for (i in 2:length(k)) {
      print(paste("k =", k[i]))
      if (fast_svd) {
        ruv3res_list[[i]] <- fastRUVIII(
          Y = t(normY),
          ctl = ctl,
          k = k[i],
          M = M,
          fullalpha = ruv3res_list[[1]]$fullalpha,
          return.info = return.info,
          propEigen = propEigen
        )
      } else {
        ruv3res_list[[i]] <- RUVIII(
          Y = t(normY),
          ctl = ctl,
          k = k[i],
          M = M,
          fullalpha = ruv3res_list[[1]]$fullalpha,
          return.info = return.info
        )
      } ## End fast_svd criterion
    } ## End length(k) == 1 criterion

    if (is.null(cell_type)) {
      cat("No cell type info, replicate matrix will be used as cell type info\n")
      cell_type <- apply(M, 1, function(x) which(x == 1))
    }

    sil_res <- do.call(cbind, lapply(ruv3res_list,
                                     FUN = function(x) {
                                       pca.data <- rsvd::rpca(x$newY, k = 10, rand = 1)
                                       # pca.data<-prcomp(x$newY)
                                       c(
                                         batch_sil(pca.data, as.numeric(as.factor(cell_type))),
                                         batch_sil(pca.data, as.numeric(as.factor(batch)), nPCs = 10)
                                       )
                                     }
    ))

    f_score <- rep(NA, ncol(sil_res))
    for (i in 1:length(k)) {
      f_score[i] <- f_measure(zeroOneScale(sil_res[1, ])[i], 1 - zeroOneScale(sil_res[2, ])[i])
    }
    names(f_score) <- k


    ruv3res <- ruv3res_list[[which.max(f_score)]]
    ruv3res$k <- which.max(f_score)
    print(sil_res)
    print(paste("optimal k:", ruv3res$k))

    plot(k, f_score, pch = 16, col = "light grey")
    lines(k, f_score)
    points(ruv3res$k, f_score[ruv3res$k], col = "red", pch = 16)
  }
  if (return_all) {
    for (i in 1:length(ruv3res_list)) {
      ruv3res_list[[i]]$newY <- t((t(ruv3res_list[[i]]$newY) * geneSdMat + geneMeanMat))
    }
    if (is.null(batch)) {
      return(ruv3res_list)
    } else {
      ruv3res_list$optimise_k <- which.max(f_score)
      return(ruv3res_list)
    }
  } else {
    ruv3res$newY <- t((t(ruv3res$newY) * geneSdMat + geneMeanMat))
    return(ruv3res)
  }
}



zeroOneScale <- function(v) {
  v <- (v+1)/2
  return(v)
}

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

f_measure <- function(celltypes, batch) {
  f <- 2 * (celltypes * batch) / (celltypes + batch)
  return(f)
}
