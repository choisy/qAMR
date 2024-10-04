#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


// Cols ----------------

// [[Rcpp::export]]
arma::rowvec Arma_colSums(const arma::mat& x) {
  return arma::sum(x, 0);
}

