#' fastRUVIII
#'
#' perform fastRUVIII
#'
#'
#' @author Yingxin Lin
#' @export


fastRUVIII <-
  function(Y, M, ctl, k = NULL, eta = NULL, propEigen = 0.1, include.intercept = TRUE,
           average = FALSE, fullalpha = NULL, return.info = FALSE, inputcheck = TRUE) {
    if (is.data.frame(Y)) {
      Y <- data.matrix(Y)
    }

    m <- nrow(Y)
    n <- ncol(Y)
    M <- ruv::replicate.matrix(M)
    ctl <- tological(ctl, n)
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
    Y <- ruv::RUV1(Y, eta, ctl, include.intercept = include.intercept)
    if (ncol(M) >= m) {
      newY <- Y
    } else if (is.null(k)) {
      ycyctinv <- solve(Y[, ctl] %*% t(Y[, ctl]))
      newY <- (M %*% solve(t(M) %*% ycyctinv %*% M) %*% (t(M) %*%
                                                           ycyctinv)) %*% Y
      fullalpha <- NULL
    }
    else if (k == 0) {
      newY <- Y
      fullalpha <- NULL
    }
    else {
      if (is.null(fullalpha)) {
        Y0 <- residop_fast(Y, M)
        ####################
        rsvd_k = min(m - ncol(M), sum(ctl), ceiling(propEigen*min(dim(Y0))))
        ####################
        fullalpha <- t(rsvd::rsvd(Y0 %*% t(Y0), k = rsvd_k, q = 2)$u[, 1:rsvd_k, drop = FALSE]) %*% Y
      }

      alpha <- fullalpha[1:min(k, nrow(fullalpha)), , drop = FALSE]
      ac <- alpha[, ctl, drop = FALSE]
      W <- Y[, ctl] %*% t(ac) %*% solve(ac %*% t(ac))
      newY <- Y - W %*% alpha
    }
    if (average) {
      newY <- ((1 / apply(M, 2, sum)) * t(M)) %*% newY
    }
    if (!return.info) {
      return(newY)
    } else {
      return(list(
        newY = newY, M = M, fullalpha = fullalpha,
        rsvd_k_options = c(
          "m-ncol(M)" = m - ncol(M),
          "sum(ctl)" = sum(ctl),
          "rsvd_prop_propEigen" = ceiling(propEigen*min(dim(Y0))))
      )
      )
    }
  }


residop_fast <-
  function(A, B) {
    return(A - B %*% solve(t(B) %*% B) %*% (t(B) %*% A))
  }

tological <-
  function(ctl, n) {
    ctl2 <- rep(FALSE, n)
    ctl2[ctl] <- TRUE
    return(ctl2)
  }


