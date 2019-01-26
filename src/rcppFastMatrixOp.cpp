// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>

//' @title fast matrix multiplication
//' @param A a matrix
//' @param B a matrix
//' @export
//' @examples
//' A = matrix(0, ncol = 1000, nrow = 1000)
//' system.time(A %*% A)
//' system.time(eigenMatMult(A, A))
// [[Rcpp::export]]
SEXP eigenMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A * B;

  return Rcpp::wrap(C);
}


//' @title fast matrix residual operator
//' @param A a matrix
//' @param B a matrix
//' @export
//' @examples
//' Y = M = diag(1, 500)
//' system.time(scMerge::eigenResidop(Y, M))
//' system.time(ruv::residop(Y, M))
// [[Rcpp::export]]
SEXP eigenResidop(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  // Eigen::MatrixXd C = A - B * (B.transpose() * B).inverse() * B.transpose() * A;
  Eigen::MatrixXd tBB_inv = (B.transpose() * B).inverse();
  Eigen::MatrixXd tBA = B.transpose() * A;
  Eigen::MatrixXd BtBB_inv = B * tBB_inv;
  Eigen::MatrixXd BtBB_inv_tBA = BtBB_inv * tBA;
  Eigen::MatrixXd C = A - BtBB_inv_tBA;
  return Rcpp::wrap(C);
}
