#include <Rcpp.h>
using namespace Rcpp;


NumericMatrix Hadamard(NumericVector v, NumericMatrix m) {
  unsigned int ncol = m.ncol();
  unsigned int nrow = m.nrow();
  NumericMatrix m2 (nrow, ncol);
  for (unsigned int i = 0; i < ncol; i++) {
    m2(_, i) = v * m(_, i);
  }
  return m2;
}


// [[Rcpp::export]]
List evolving1(double dx, double alpha, double mu, NumericMatrix J,
               NumericVector ps, NumericVector dt_val, NumericVector ks,
               NumericVector mr_bact_pop_in) {
  int nt = dt_val.size() + 1;
  int nx = ps.size();
  double mr_pop_size;
  NumericMatrix bact_pop (nx, nt);
  NumericVector pop_size (nt);
  NumericVector mr_bact_pop (nx);
  mr_bact_pop = clone(mr_bact_pop_in);
  bact_pop(_, 0) = mr_bact_pop;
  mr_pop_size = sum(mr_bact_pop);
  pop_size(0) = mr_pop_size;
  for(int i = 1; i < nt; i++) {
    mr_bact_pop = (mr_bact_pop + dt_val(i - 1) * (colSums(Hadamard(ps * mr_bact_pop, J)) *
      dx / pow(1.0 + mr_pop_size, alpha))) / (1.0 + dt_val(i - 1) * (mu + ks));
    bact_pop(_, i) = mr_bact_pop;
    mr_pop_size = sum(mr_bact_pop);
    pop_size(i) = mr_pop_size;
  }
  return List::create(bact_pop, pop_size);
}


// [[Rcpp::export]]
List evolving2(double dx, double alpha, double mu,
               NumericVector ps, NumericVector dt_val, NumericVector ks,
               NumericVector mr_bact_pop_in) {
  int nt = dt_val.size() + 1;
  int nx = ps.size();
  double mr_pop_size;
  NumericMatrix bact_pop (nx, nt);
  NumericVector pop_size (nt);
  NumericVector mr_bact_pop (nx);
  mr_bact_pop = clone(mr_bact_pop_in);
  bact_pop(_, 0) = mr_bact_pop;
  mr_pop_size = sum(mr_bact_pop);
  pop_size(0) = mr_pop_size;
  for(int i = 1; i < nt; i++) {
    mr_bact_pop = (mr_bact_pop + dt_val(i - 1) *
      (ps * mr_bact_pop / pow(1.0 + mr_pop_size, alpha))) / (1.0 + dt_val(i - 1) * (mu + ks));
    bact_pop(_, i) = mr_bact_pop;
    mr_pop_size = sum(mr_bact_pop);
    pop_size(i) = mr_pop_size;
  }
  return List::create(bact_pop, pop_size);
}