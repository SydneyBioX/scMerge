#' @title Simulate a simple matrix to test the speed of scMerge
#' @param m Number of observations
#' @param n Number of features
#' @param nc Number of negative controls
#' @param nRep Number of replicates
#' @param lambda Rate parameter for random Poisson generation
#' @param sce If \code{TRUE}, returns a SingleCellExperiment object
#' @description This function is designed to generate Poisson-random-variable data matrix
#' to test on the speed and performance of the internal algorithms of scMerge.
#' It is not intended to be used by end0-users.
#' @export
#' @examples
#' L = ruvSimulate(m = 500, n = 20000, nRep = 10, lambda = 0.1)
#' Y = L$Y; M = L$M;
#'
#' system.time(scMerge::eigenResidop(Y, M))
#' system.time(ruv::residop(Y, M))
#'
#'
#'

ruvSimulate = function(m = 100, n = 1e5, nc = 1e3, nRep = 50, lambda = 0.1, sce = FALSE){
  # m Number of observations
  # n Number of features
  # nc Number of negative controls
  p = 1
  k = 20
  ctl = rep(FALSE, n)
  ctl[1:nc] = TRUE
  X = matrix(c(rep(0,floor(m/2)), rep(1,ceiling(m/2))), m, p) ## Design matrix is simple two conditions
  beta = matrix(stats::rpois(p*n, lambda = lambda), p, n)
  beta[,ctl] = 0
  W = matrix(stats::rpois(m*k, lambda = lambda),m,k)
  alpha = matrix(stats::rpois(k*n, lambda = lambda),k,n)
  epsilon = matrix(stats::rpois(m*n, lambda = lambda),m,n)
  Y = X %*% beta + W %*% alpha + epsilon

  # Define patientID and sampleDate
  # patientID = paste("patient", rep_len(1:nRep, length.out = m), sep = "")
  # #print(patientID)
  # sampleDate = paste("Date", rep_len(c(1:4), length.out = m), sep = "")
  # batch = paste(patientID, sampleDate, sep = "_")
  # M = ruv::replicate.matrix(data.frame(patientID, sampleDate))

  cellType = paste("cellType", rep_len(1:nRep, length.out = m), sep = "")
  dataSource = paste("dataSource", rep_len(c(1:4), length.out = m), sep = "")
  # batch = paste(cellType, dataSource, sep = "_")
  M = ruv::replicate.matrix(data.frame(cellType, dataSource))

  if(sce){
    result = SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = t(Y),
                    logcounts = log2(t(Y) + 1L)),
      colData = data.frame(cellType, dataSource)
    )
  } else {
    result = list(Y = Y,
                  ctl = ctl,
                  M = M,
                  cellType = cellType,
                  dataSource = dataSource)
  }
  return(result)
}
