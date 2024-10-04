#include <Rcpp.h>
using namespace Rcpp;

NumericVector colSums(NumericMatrix x) {
  int nr = x.nrow(), nc = x.ncol();
  NumericVector ans(nc);
  for (int j = 0; j < nc; j++) {
    double sum = 0.0;
    for (int i = 0; i < nr; i++) {
      sum += x(i, j);
    }
    ans[j] = sum;
  }
  return ans;
}


// [[Rcpp::export]]
List evolving1(double dx, double alpha, double mu,
               NumericVector J, NumericVector ps, NumericVector dt_val, NumericVector ks,
               NumericVector mr_bact_pop) {
  int nt = dt_val.size();
  int nx = ps.size();
  double mr_pop_size;
  NumericMatrix bact_pop (nx, nt);
  NumericVector pop_size (nt);
  bact_pop(_, 0) = mr_bact_pop;
  pop_size(0) = sum(mr_bact_pop);
  for(int i = 0; i < nt; i++) {
    mr_bact_pop = (mr_bact_pop + dt_val(i) * colSums(J *  (ps % bact_pop)) * dx /
      (1 + pop_size(i))^alpha) / (1 + dt_val(i) * (mu + ks));
    bact_pop(_, i + 1) = mr_bact_pop;
    pop_size(i + 1) = sum(mr_bact_pop);
  }
  return List::create(bact_pop, pop_size);
}


