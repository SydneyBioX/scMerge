#' @title Simulate a simple matrix
#' @param m Number of observations
#' @param n Number of features
#' @param nc Number of negative controls
#' @param nRep Number of replicates
#' @export
#' @examples
#' L = ruvSimulate(m = 700, n = 10000, nRep = 10, lambda = 0.1)
#' Y = L$Y; M = L$M; ctl = L$ctl; batch = L$batch
#'
#' microbenchmark::microbenchmark(
#' improved1 = scMerge::fastRUVIII(Y = Y, M = M, ctl = ctl, k = 20, fast_svd = FALSE),
#' old = ruv::RUVIII(Y = Y, M = M, ctl = ctl, k = 20),
#' times = 1)
#'
#' microbenchmark::microbenchmark(scMerge::scRUVIII(Y = Y, M = M, ctl = ctl, k = 20, batch = batch), times = 1)
#' microbenchmark::microbenchmark(fakePackage::scRUVIII(Y = Y, M = M, ctl = ctl, k = 20, batch = batch), times = 1)
#'
ruvSimulate = function(m = 100, n = 1e5, nc = 1e3, nRep = 50, lambda = 0.1){
  # m Number of observations
  # n Number of features
  # nc Number of negative controls
  p = 1
  k = 20
  ctl = rep(FALSE, n)
  ctl[1:nc] = TRUE
  X = matrix(c(rep(0,floor(m/2)), rep(1,ceiling(m/2))), m, p)
  beta = matrix(rpois(p*n, lambda = lambda), p, n)
  beta[,ctl] = 0
  W = matrix(rpois(m*k, lambda = lambda),m,k)
  alpha = matrix(rpois(k*n, lambda = lambda),k,n)
  epsilon = matrix(rpois(m*n, lambda = lambda),m,n)
  Y = X%*%beta + W%*%alpha + epsilon

  # Define patientID and sampleDate
  patientID = paste("patient", rep_len(1:nRep, length.out = m), sep="")
  #print(patientID)
  sampleDate = paste("Date", rep_len(c(1:4), length.out = m), sep="")
  batch = paste(patientID, sampleDate, sep = "_")
  M = ruv::replicate.matrix(data.frame(patientID, sampleDate))

  result = list(Y = Y,
                ctl = ctl,
                M = M,
                batch = batch)
  return(result)
}
