#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat kron_function(arma::mat A, arma::mat B) {
  return kron(A, B);
}
