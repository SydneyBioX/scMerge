// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>

//' @title fast matrix multiplication
//' @param A a matrix
//' @param B a matrix
//' @export
// [[Rcpp::export]]
SEXP eigenMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A * B;

  return Rcpp::wrap(C);
}


//' @title fast matrix inverse
//' @param A a matrix
//' @export
// [[Rcpp::export]]
SEXP eigenMatInverse(const Eigen::Map<Eigen::MatrixXd> A){
  Eigen::MatrixXd Ainv = A.inverse();

  return Rcpp::wrap(Ainv);
}


//' @title fast matrix inverse
//' @param A a matrix
//' @param B a matrix
//' @export
// [[Rcpp::export]]
SEXP eigenResidop2(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  // Eigen::MatrixXd C = A - B * (B.transpose() * B).inverse() * B.transpose() * A;
  Eigen::MatrixXd tBB_inv = (B.transpose() * B).inverse();
  Eigen::MatrixXd tBA = B.transpose() * A;
  Eigen::MatrixXd BtBB_inv = B * tBB_inv;
  Eigen::MatrixXd BtBB_inv_tBA = BtBB_inv * tBA;
  Eigen::MatrixXd C = A - BtBB_inv_tBA;
  return Rcpp::wrap(C);
}
