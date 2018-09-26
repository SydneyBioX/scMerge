// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>


// [[Rcpp::export]]
SEXP eigenMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A * B;
  
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenMatInverse(const Eigen::Map<Eigen::MatrixXd> A){
  Eigen::MatrixXd Ainv = A.inverse();
  
  return Rcpp::wrap(Ainv);
}