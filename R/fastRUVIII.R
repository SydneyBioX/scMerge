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
#' @param BPPARAM A \code{BiocParallelParam} class object from the \code{BiocParallel} package is used. Default is SerialParam().
#' @param BSPARAM A \code{BiocSingularParam} class object from the \code{BiocSingular} package is used. Default is ExactParam().
#' @param svd_k If BSPARAM is set to \code{RandomParam} or \code{IrlbaParam} class from \code{BiocSingular} package, then 
#' \code{svd_k} will be used to used to reduce the computational cost of singular value decomposition. Default to 50.
#' @param include.intercept When eta is specified (not NULL) but does not already include an intercept term, this will automatically include one.
#' See ruv::RUVIII for details.
#' @param average Average replicates after adjustment. See ruv::RUVIII for details.
#' @param fullalpha Not used. Please ignore. See ruv::RUVIII for details.
#' @param return.info Additional information relating to the computation of normalised matrix. We recommend setting this to true.
#' @param inputcheck We recommend setting this to true.
#' @useDynLib scMerge, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom DelayedArray t
#' @importFrom BiocSingular runExactSVD
#' @export
#' @return
#' A normalised matrix of the same dimensions as the input matrix Y.
#' @examples
#' L = ruvSimulate(m = 200, n = 500, nc = 400, nCelltypes = 3, nBatch = 2, lambda = 0.1, sce = FALSE)
#' Y = L$Y; M = L$M; ctl = L$ctl
#' improved1 = scMerge::fastRUVIII(Y = Y, M = M, ctl = ctl, 
#' k = 20, BSPARAM = BiocSingular::ExactParam())
#' improved2 = scMerge::fastRUVIII(Y = Y, M = M, ctl = ctl, 
#' k = 20, BSPARAM = BiocSingular::RandomParam(), svd_k = 50)
#' old = ruv::RUVIII(Y = Y, M = M, ctl = ctl, k = 20)
#' all.equal(improved1, old)
#' all.equal(improved2, old)


fastRUVIII <- function(Y, M, ctl, k = NULL, eta = NULL,
    svd_k = 50, include.intercept = TRUE, average = FALSE, 
    BPPARAM = SerialParam(), BSPARAM = ExactParam(),
    fullalpha = NULL, return.info = FALSE, inputcheck = TRUE) {
    
    m <- nrow(Y)
    n <- ncol(Y)
    M <- ruv::replicate.matrix(M)
    ctl <- tological(ctl, n)
    
    ## Check the inputs
    if (inputcheck) {
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
    
    if (class(BSPARAM) != "ExactParam") {
        svd_k <- min(m - ncol(M), sum(ctl), svd_k, na.rm = TRUE)
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
                
                if(class(Y) == "matrix"){
                    Y0 <- eigenResidop(Y, M)
                } else if (class(Y) == "dgeMatrix"){
                    Y0 <- eigenResidop(as.matrix(Y), M)
                } else {
                    Y0 <- ruv::residop(Y, M)
                }
            
                svdObj <- BiocSingular::runSVD(x = Y0, k = svd_k, BPPARAM = BPPARAM, BSPARAM = BSPARAM)
                
                fullalpha <- t(svdObj$u[, seq_len(svd_k), drop = FALSE]) %*% Y
            }  ## End is.null(fullalpha)
        ###############
        alpha <- fullalpha[seq_len(min(k, nrow(fullalpha))), 
            , drop = FALSE]
        ac <- alpha[, ctl, drop = FALSE]
        W <- Y[, ctl] %*% t(ac) %*% solve(ac %*% t(ac))
        newY <- Y - W %*% alpha
    }  ## End else(ncol(M) >= m | k == 0)
    
    ## If the users want to get all the informations relating to
    ## the RUV, it can be done here.
    if (!return.info) {
        return(newY)
    } else {
        return(list(newY = newY, M = M, fullalpha = fullalpha))
    }
}
############################ 
tological <- function(ctl, n) {
    ctl2 <- rep(FALSE, n)
    ctl2[ctl] <- TRUE
    return(ctl2)
}
