#' scRUVIII
#'
#' perform scRUVIII
#'
#'
#' @author Yingxin Lin
#' @export


# A function to perform loacation/scale adjustement to data as the input of RUVIII
# which also provides the option to select optimal RUVk according to the silhouette coefficient
scRUVIII <- function(Y = Y, M = M, ctl = ctl,
                             fullalpha = NULL,
                             k = k, return.info = T,
                             cell_type = NULL,
                             batch = NULL,
                             return_all = F,
                             fast_svd = F) {
  Y <- t(Y)
  # geneMeanMat<-matrix(rep(rowMeans(Y),ncol(Y)),ncol=ncol(Y))
  # geneSdMat<-matrix(rep(apply(Y-geneMeanMat,1,sd),ncol(Y)),ncol=ncol(Y))
  # normY<-(Y-geneMeanMat)/geneSdMat

  if (is.null(batch)) {
    warning("No batch info!")
    return(NULL)
  }

  scale_res <- standardize(Y, batch)
  normY <- scale_res$s.data
  geneSdMat <- sqrt(scale_res$stand.var) %*% t(rep(1, ncol(Y)))
  geneMeanMat <- scale_res$stand.mean

  if (length(k) == 1) {
    if (fast_svd) {
      ruv3res <- fastRUVIII(
        Y = t(normY),
        ctl = ctl,
        k = k,
        M = M,
        fullalpha = fullalpha,
        return.info = return.info
      )
    } else {
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
  } else {
    ruv3res_list <- list()
    print(paste("k =", k[1]))
    if (fast_svd) {
      ruv3res_list[[1]] <- fastRUVIII(
        Y = t(normY),
        ctl = ctl,
        k = k[1],
        M = M,
        fullalpha = fullalpha,
        return.info = return.info
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
    }

    for (i in 2:length(k)) {
      print(paste("k =", k[i]))
      if (fast_svd) {
        ruv3res_list[[i]] <- fastRUVIII(
          Y = t(normY),
          ctl = ctl,
          k = k[i],
          M = M,
          fullalpha = ruv3res_list[[1]]$fullalpha,
          return.info = return.info
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
      }
    }

    if (is.null(cell_type)) {
      cat("No cell type info, replicate matrix will be used as cell type info\n")
      cell_type <- apply(M, 1, function(x) which(x == 1))
    }

    sil_res <- do.call(cbind, lapply(ruv3res_list,
                                     FUN = function(x) {
                                       pca.data <- rpca(x$newY, k = 50, rand = 1)
                                       # pca.data<-prcomp(x$newY)
                                       c(
                                         batch_sil(pca.data, as.numeric(as.factor(cell_type))),
                                         batch_sil(pca.data, as.numeric(as.factor(batch)), nPCs = 10)
                                       )
                                     }
    ))
    # score<-minMaxScale(k*sil_res[1,])-minMaxScale(rev(k)*sil_res[2,])
    # f_score<-minMaxScale(sil_res[1,])-minMaxScale(sil_res[2,])
    f_score <- rep(NA, ncol(sil_res))
    for (i in 1:length(k)) {
      f_score[i] <- f_measure(minMaxScale(sil_res[1, ])[i], 1 - minMaxScale(sil_res[2, ])[i])
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



minMaxScale <- function(v) {
  v <- (v - min(v)) / (max(v) - min(v))
  return(v)
}

standardize <- function(exprsMat, batch) {
  num_cell <- ncol(exprsMat)
  batch <- as.factor(batch)
  grand.mean <- matrix(rowMeans(exprsMat), nrow = 1)
  stand.mean <- t(grand.mean) %*% t(rep(1, num_cell))
  design <- model.matrix(~-1 + batch)
  B.hat <- solve(crossprod(design), tcrossprod(t(design), as.matrix(exprsMat)))
  var.pooled <- ((exprsMat - t(design %*% B.hat))^2) %*% rep(1 / num_cell, num_cell)
  s.data <- (exprsMat - stand.mean) / (sqrt(var.pooled) %*% t(rep(1, num_cell)))
  return(res = list(s.data = s.data, stand.mean = stand.mean, stand.var = var.pooled))
}


