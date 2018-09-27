#' @title Simulate a simple matrix
#' @param m Number of observations
#' @param n Number of features
#' @param nc Number of negative controls
#' @param nRep Number of replicates
#' @export
#' 
ruvSimulate = function(m = 100, n = 1e5, nc = 1e3, nRep = 50){
  # m Number of observations
  # n Number of features
  # nc Number of negative controls
  p = 1
  k = 20
  ctl = rep(FALSE, n)
  ctl[1:nc] = TRUE
  X = matrix(c(rep(0,floor(m/2)), rep(1,ceiling(m/2))), m, p)
  beta = matrix(rpois(p*n, 1), p, n)
  beta[,ctl] = 0
  W = matrix(rpois(m*k, 1),m,k)
  alpha = matrix(rpois(k*n,1),k,n)
  epsilon = matrix(rpois(m*n,1),m,n)
  Y = X%*%beta + W%*%alpha + epsilon
  
  # Define patientID and sampleDate
  patientID = paste("patient", rep_len(1:nRep, length.out = m), sep="")
  #print(patientID)
  sampleDate = paste("Date", rep_len(c(1:4), length.out = m), sep="")
  
  M = ruv::replicate.matrix(data.frame(patientID, sampleDate))
  
  result = list(Y = Y, 
                ctl = ctl,
                M = M)
  return(result)
}