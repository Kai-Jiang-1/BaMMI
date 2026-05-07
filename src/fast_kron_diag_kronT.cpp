#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat fast_kron_diag_kronT(const arma::mat& X,
                               const arma::mat& A,
                               const arma::vec& w) {
  int n = X.n_rows;
  int p = X.n_cols;
  int q = A.n_rows;

  arma::mat M(n * q, n * q, arma::fill::zeros);

  int idx = 0;
  for (int i = 0; i < p; ++i) {
    arma::vec xi  = X.col(i);
    arma::mat xxt = xi * xi.t();

    for (int j = 0; j < q; ++j) {
      arma::vec aj  = A.col(j);
      arma::mat aat = aj * aj.t();

      M += w(idx) * arma::kron(xxt, aat);
      ++idx;
    }
  }

  return M;
}
