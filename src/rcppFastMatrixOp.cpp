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
