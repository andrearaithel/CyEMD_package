#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
double emdC(Rcpp::NumericVector a, Rcpp::NumericVector b) {
  int n = a.size();
  Rcpp::NumericVector dist = Rcpp::NumericVector(n);
  double emd = 0;
  for(int i = 0; i < (n - 1); ++i) {
    dist[i + 1] = a[i] - b[i] + dist[i];
  }
  dist = abs(dist);
  for (auto& d : dist)
    emd += d;
  return emd;
}
