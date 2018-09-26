#' fastRUVIII
#'
#' perform fastRUVIII
#' @author Yingxin Lin, John Ormerod, Kevin Wang
#' @param Y The unnormalised SC data. A m by n matrix, where m is the number of observations and n is the number of features.
#' @param M The replicate mapping matrix.
#' The mapping matrix has m rows (one for each observation), and each column represents a set of replicates. The (i, j)-th entry of the mapping matrix is 1 if the i-th observation is in replicate set j, and 0 otherwise.
#' See ruv::RUVIII for more details.
#' @param ctl An index vector to specify the negative controls. Either a logical vector of length n or a vector of integers.
#' @param k The number of unwanted factors to remove. This is inherited from the ruvK argument from the scMerge::scMerge function.
#' @param eta Gene-wise (as opposed to sample-wise) covariates. See ruv::RUVIII for details.
#' @param rsvd_prop rsvd_prop is used to control for the accuracy of the randomized SVD computation.
#' If a lower value is used on a lower dimensional data (say < 1000 cell) will potentially yield a less accurate computed result but with a gain in speed.
#' The default of 0.1 tends to achieve a balance between speed and accuracy.
#' @param include.intercept When eta is specified (not NULL) but does not already include an intercept term, this will automatically include one.
#' See ruv::RUVIII for details.
#' @param average Average replicates after adjustment. See ruv::RUVIII for details.
#' @param fullalpha Not used. Please ignore. See ruv::RUVIII for details.
#' @param return.info Additional information relating to the computation of normalised matrix. We recommend setting this to true.
#' @param inputcheck We recommend setting this to true.
#' @export


fastRUVIII <-
  function(Y, M, ctl, k = NULL, eta = NULL, rsvd_prop = 0.1, include.intercept = TRUE,
           average = FALSE, fullalpha = NULL, return.info = FALSE, inputcheck = TRUE) {

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
        warning("m is greater than n!  This is not a problem itself, but may indicate that you need to transpose your data matrix.  Please ensure that rows correspond to observations (e.g. microarrays) and columns correspond to features (e.g. genes).")
      }
      if (sum(is.na(Y)) > 0) {
        warning("Y contains missing values.  This is not supported.")
      }
      if (sum(Y == Inf, na.rm = TRUE) + sum(Y == -Inf, na.rm = TRUE) >
          0) {
        warning("Y contains infinities.  This is not supported.")
      }
    }


    ## RUV1 is a reprocessing step for RUVIII
    Y <- ruv::RUV1(Y, eta, ctl, include.intercept = include.intercept)


    ## m represent the number of samples/observations
    ## ncol(M) represent the number of replicates
    ## If the replicate matrix is such that we have more replicates than samples, then RUV3 is not appropriate, thus, we return the
    ## Original input matrix
    if (ncol(M) >= m) {
      newY <- Y
    } else if (is.null(k)) { ## If there isn't a value for RUVk, we will compute an alternative type of normalisation. In scMerge, RUVk is a forced input.
      ycyctinv <- solve(Y[, ctl] %*% t(Y[, ctl]))
      newY <- (M %*% solve(t(M) %*% ycyctinv %*% M) %*% (t(M) %*%
                                                           ycyctinv)) %*% Y
      fullalpha <- NULL
    }
    else if (k == 0) { ## If RUVk = 0, there is no need for normalisation.
      newY <- Y
      fullalpha <- NULL
    }
    else {
      if (is.null(fullalpha)) { ## The main RUVIII process
        ## Applies the residual operator of a matrix M to a matrix Y
        Y0 <- eigenResidop(Y, M)
        ## In RSVD, we need to specify a paritial decomposition parameter, rsvd_k.
        ## In the usual RUVIII procedure, we can set this to min(m - ncol(M), sum(ctl))
        ## But it was found that we can lower this parameter to ceiling(rsvd_prop*m), if we want to speed up computation.
        rsvd_k = min(m - ncol(M), sum(ctl), ceiling(rsvd_prop*m), na.rm = TRUE)
        ####################
        svdObj = rsvd::rsvd(eigenMatMult(Y0, t(Y0)), k = rsvd_k)
        fullalpha = eigenMatMult(t(svdObj$u[, 1:rsvd_k, drop = FALSE]), Y)
      }

      ## alpha matrix is a subset of of the fullalpha matrix
      alpha <- fullalpha[1:min(k, nrow(fullalpha)), , drop = FALSE]
      ac <- alpha[, ctl, drop = FALSE]
      W <- Y[, ctl] %*% t(ac) %*% solve(ac %*% t(ac))
      newY <- Y - W %*% alpha
    }
    if (average) { ## If average over the replicates is needed. This is ignored in scMerge.
      newY <- ((1 / apply(M, 2, sum)) * t(M)) %*% newY
    }

    ## If the users want to get all the informations relating to the RUV, it can be done here.
    if (!return.info) {
      return(newY)
    } else {
      return(list(
        newY = newY, M = M, fullalpha = fullalpha,
        rsvd_k_options = c(
          "m-ncol(M)" = m - ncol(M),
          "sum(ctl)" = sum(ctl),
          "rsvd_prop_rsvd_prop" = ceiling(rsvd_prop*min(dim(Y0))))
      )
      )
    }
  }

############################
residop_fast <-
  function(A, B) {
    return(A - B %*% (solve(t(B) %*% B) %*% (t(B) %*% A)))
  }
############################
tological <-
  function(ctl, n) {
    ctl2 <- rep(FALSE, n)
    ctl2[ctl] <- TRUE
    return(ctl2)
  }
