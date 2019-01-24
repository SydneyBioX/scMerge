#' @title RUVg function for single cell
#' @description Modified based on RUV2 from package ruv and RUVg from package RUVseq function
#' @import ruv
#' @author Yingxin Lin

scRUVg <-function (Y,  ctl, k, Z = 1, eta = NULL, include.intercept = TRUE,
                   fullW = NULL, svdyc = NULL) {
  if (is.data.frame(Y)){
    Y = data.matrix(Y)
  }

  m = nrow(Y)
  n = ncol(Y)
  if (is.numeric(Z)) {
    if (length(Z) == 1) {
      Z = matrix(1, m, 1)
    }
  }

  if (!is.null(Z)) {
    Z = ruv::design.matrix(Z, name = "Z", include.intercept = include.intercept)
    q = ncol(Z)
  }else{
    q = 0
  }

  ctl2 <- rep(FALSE, n)
  ctl2[ctl] <- TRUE
  ctl <- ctl2

  if (k > sum(ctl)) {
    stop("k must not be larger than the number of controls")
  }

  Y0 = ruv::RUV1(Y, eta, ctl, include.intercept = include.intercept)

  if (q > 0) {
    Y0 = residop(Y0, Z)
  }
  if (is.null(fullW)) {
    if (is.null(svdyc)){
      svdyc = svd(Y0[, ctl, drop = FALSE]%*%t(Y0[, ctl, drop = FALSE]))
      fullW = svdyc$u[, seq_len(min((m - q), sum(ctl))), drop = FALSE]
    }
  }

  W = alpha =  NULL

  W = fullW[, seq_len(k), drop = FALSE]
  alpha = solve(t(W) %*% W) %*% t(W) %*% Y0

  newY = Y - W %*% alpha

  return(list(newY = newY, W = W, alpha = alpha))
}


