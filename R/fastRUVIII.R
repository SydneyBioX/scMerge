#' @title A fast version of the ruv::RUVIII algorithm
#' @description Perform a fast version of the ruv::RUVIII algorithm for scRNA-Seq data noise estimation
#'
#' @author Yingxin Lin, John Ormerod, Kevin Wang
#' @param Y The unnormalised scRNA-Seq data matrix. 
#' A m by n matrix, where m is the number of observations and n is the number of features.
#' @param M The replicate mapping matrix.
#' The mapping matrix has m rows (one for each observation), and each column represents a set of replicates. 
#' The (i, j)-th entry of the mapping matrix is 1 if the i-th observation is in replicate set j, and 0 otherwise.
#' See ruv::RUVIII for more details.
#' @param ctl An index vector to specify the negative controls. Either a logical vector of length n or a vector of integers.
#' @param k The number of unwanted factors to remove. This is inherited from the ruvK argument from the scMerge::scMerge function.
#' @param eta Gene-wise (as opposed to sample-wise) covariates. See ruv::RUVIII for details.
#' @param fast_svd If \code{TRUE}, fast algorithms will be used for singular value decomposition calculation via the \code{irlba} and \code{rsvd} packages. 
#' We recommend using this option when the number of cells is large (e.g. more than 1000 cells).
#' @param rsvd_prop If \code{fast_svd = TRUE}, then \code{rsvd_prop} will be used to used to reduce the computational cost of randomised singular value decomposition. 
#' We recommend setting this number to less than 0.25 to achieve a balance between numerical accuracy and computational costs.
#' @param include.intercept When eta is specified (not NULL) but does not already include an intercept term, this will automatically include one.
#' See ruv::RUVIII for details.
#' @param average Average replicates after adjustment. See ruv::RUVIII for details.
#' @param fullalpha Not used. Please ignore. See ruv::RUVIII for details.
#' @param return.info Additional information relating to the computation of normalised matrix. We recommend setting this to true.
#' @param inputcheck We recommend setting this to true.
#' @useDynLib scMerge, .registration = TRUE
#' @importFrom rsvd rsvd 
#' @importFrom Rcpp sourceCpp
#' @import RcppEigen
#' @export
#' @return
#' A normalised matrix of the same dimensions as the input matrix Y.
#' @examples
#' L = ruvSimulate(m = 200, n = 500, nc = 400, nCelltypes = 3, nBatch = 2, lambda = 0.1, sce = FALSE)
#' Y = L$Y; M = L$M; ctl = L$ctl
#' improved1 = scMerge::fastRUVIII(Y = Y, M = M, ctl = ctl, k = 20, fast_svd = FALSE)
#' improved2 = scMerge::fastRUVIII(Y = Y, M = M, ctl = ctl, k = 20, fast_svd = TRUE, rsvd_prop = 0.1)
#' old = ruv::RUVIII(Y = Y, M = M, ctl = ctl, k = 20)
#' all.equal(improved1, old)
#' all.equal(improved2, old)


fastRUVIII <- function(Y, M, ctl, k = NULL, eta = NULL, fast_svd = FALSE, 
    rsvd_prop = 0.1, include.intercept = TRUE, average = FALSE, 
    fullalpha = NULL, return.info = FALSE, inputcheck = TRUE) {
    
    if (is.data.frame(Y)) {
        Y <- data.matrix(Y)
    }
    
    m <- nrow(Y)
    n <- ncol(Y)
    M <- ruv::replicate.matrix(M)
    ctl <- tological(ctl, n)
    
    ## Check the inputs
    if (inputcheck) {
        if (m > n) {
            message("m is greater than n! 
            This is not a problem itself, but may indicate that you need to transpose your data matrix. 
                    Please ensure that rows correspond to observations (e.g. microarrays) and columns correspond to features (e.g. genes).")
        }
        if (sum(is.na(Y)) > 0) {
            stop("Y contains missing values.  This is not supported.")
        }
        if (sum(Y == Inf, na.rm = TRUE) + sum(Y == -Inf, na.rm = TRUE) > 
            0) {
            stop("Y contains infinities.  This is not supported.")
        }
    }
    
    
    ## RUV1 is a reprocessing step for RUVIII
    Y <- ruv::RUV1(Y, eta, ctl, include.intercept = include.intercept)
    
    if (fast_svd) {
        svd_k <- min(m - ncol(M), sum(ctl), ceiling(rsvd_prop * 
            m), na.rm = TRUE)
    } else {
        svd_k <- min(m - ncol(M), sum(ctl), na.rm = TRUE)
    }
    
    
    ## m represent the number of samples/observations ncol(M)
    ## represent the number of replicates If the replicate matrix
    ## is such that we have more replicates than samples, then
    ## RUV3 is not appropriate, thus, we return the Original input
    ## matrix
    if (ncol(M) >= m | k == 0) {
        newY <- Y
        fullalpha <- NULL
    } else {
        
        if (is.null(fullalpha)) 
            {
                ## The main RUVIII process Applies the residual operator of a
                ## matrix M to a matrix Y Y0 has the same dimensions as Y,
                ## i.e. m rows (observations) and n columns (genes).
                
                Y0 <- eigenResidop(Y, M)
                matToDecomp <- Y0
                
                if (max(dim(Y0))/min(dim(Y0)) >= 5) {
                  matToDecomp <- eigenMatMult(Y0, t(Y0))
                }
                
                
                if (fast_svd) {
                  svdObj <- rsvd::rsvd(matToDecomp, k = svd_k)
                } else {
                  svdObj <- base::svd(matToDecomp)
                }  ## End if(fast_svd)
                
                # fullalpha <- eigenMatMult(t(svdObj$u[, 1:svd_k, drop =
                # FALSE]), Y)
                fullalpha <- eigenMatMult(t(svdObj$u[, seq_len(svd_k), 
                  drop = FALSE]), Y)
            }  ## End is.null(fullalpha)
        ###############
        alpha <- fullalpha[seq_len(min(k, nrow(fullalpha))), 
            , drop = FALSE]
        ac <- alpha[, ctl, drop = FALSE]
        W <- Y[, ctl] %*% t(ac) %*% solve(ac %*% t(ac))
        newY <- Y - eigenMatMult(W, alpha)
    }  ## End else(ncol(M) >= m | k == 0)
    
    if (average) {
        ## If average over the replicates is needed. This is ignored
        ## in scMerge.
        newY <- ((1/apply(M, 2, sum)) * t(M)) %*% newY
    }
    
    ## If the users want to get all the informations relating to
    ## the RUV, it can be done here.
    if (!return.info) {
        return(newY)
    } else {
        return(list(newY = newY, M = M, fullalpha = fullalpha, 
            rsvd_k_options = c(`m-ncol(M)` = m - ncol(M), `sum(ctl)` = sum(ctl), 
                svd_k = svd_k)))
    }
}
############################ 
tological <- function(ctl, n) {
    ctl2 <- rep(FALSE, n)
    ctl2[ctl] <- TRUE
    return(ctl2)
}
