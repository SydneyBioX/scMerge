#' @title fast residop
#' @import Rcpp
#' @importFrom Rcpp sourceCpp
#' @import RcppEigen
#' @useDynLib scMerge
#' @export

eigenResidop = function(A, B){

  tBB = eigenMatMult(t(B), B)
  tBA = eigenMatMult(t(B), A)
  tBB_inv = eigenMatInverse(tBB)
  BtBB_inv = eigenMatMult(B, tBB_inv)
  BtBB_inv_tBA = eigenMatMult(BtBB_inv, tBA)

  return(A - BtBB_inv_tBA)
}
