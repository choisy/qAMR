#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix vectmat(NumericMatrix m, NumericVector v) {
  return m % v;
}


