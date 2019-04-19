#' @title RUVg function for single cell
#' @description Modified based on RUV2 from package ruv and RUVg from package RUVseq function
#' (see these function's documentations for full documentations and usage)
#' @param Y The data. A m by n matrix, where m is the number of observations and n is the number of features.
#' @param ctl index vector to specify the negative controls.
#' @param k The number of unwanted factors to use.
#' @param Z Any additional covariates to include in the model.
#' @param eta Gene-wise (as opposed to sample-wise) covariates.
#' @param include.intercept Applies to both Z and eta. When Z or eta (or both) is
#' specified (not NULL) but does not already include an intercept term, this will automatically include one.
#' If only one of Z or eta should include an intercept, this variable should be set to FALSE,
#' and the intercept term should be included manually where desired.
#' @param fullW Can be included to speed up execution. Is returned by previous calls of scRUVg
#' @param svdyc Can be included to speed up execution. For internal use; please use fullW instead.
#' @import ruv
#' @export
#' @author Yingxin Lin, Kevin Wang
#' @return A list consists of:
#' \itemize{
#' \item A matrix newY, the normalised matrix,
#' \item A matrix W, the unwanted variation matrix, and ;
#' \item A matrix alpha, this corresponding coefficient matrix for W.
#' }
#' @examples
#' L = scMerge::ruvSimulate(m = 800, n = 1000, nc = 50, nCelltypes = 10)
#' Y = L$Y; M = L$M; ctl = L$ctl
#' ruvgRes = scMerge::scRUVg(Y = Y, ctl = ctl, k = 20)


scRUVg <- function(Y, ctl, k, Z = 1, eta = NULL, include.intercept = TRUE, fullW = NULL, 
    svdyc = NULL) {
    if (is.data.frame(Y)) {
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
    } else {
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
        Y0 = eigenResidop(Y0, Z)
    }
    if (is.null(fullW)) {
        if (is.null(svdyc)) {
            Y0ctl = Y0[, ctl, drop = FALSE]
            matToDecomp = Y0ctl
            if (max(dim(matToDecomp))/min(dim(matToDecomp)) >= 5) {
                matToDecomp <- eigenMatMult(Y0ctl, t(Y0ctl))
            }
            svdyc = svd(matToDecomp)
            fullW = svdyc$u[, seq_len(min((m - q), sum(ctl))), drop = FALSE]
        }
    }
    
    W = alpha = NULL
    
    W = fullW[, seq_len(k), drop = FALSE]
    alpha = solve(t(W) %*% W) %*% t(W) %*% Y0
    
    newY = Y - W %*% alpha
    
    return(list(newY = newY, W = W, alpha = alpha))
}


