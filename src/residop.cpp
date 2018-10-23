#include <Rcpp.h>
// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;
using namespace Rcpp;

//' residop
//'
//' @usage residop(A, B)
//' @param A The first matrix
//' @param B The second matrix
//' @return The result
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd residop(Eigen::MatrixXd A, Eigen::MatrixXd B)
{
    return A - B * (B.transpose() * B).inverse() * B.transpose() * A;
}
