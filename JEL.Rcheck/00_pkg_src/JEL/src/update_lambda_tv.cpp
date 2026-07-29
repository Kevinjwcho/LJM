#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List update_lambda_tv(List Es_exp, NumericVector l0, int n) {
  List lambda(n);
  
  for (int i = 0; i < n; i++) {
    NumericVector y = Es_exp[i];
    int len_t = y.size();
    // std::fill(y.begin(), y.end(), Es_exp[i]);
    
    if (len_t == l0.size()) {
      lambda[i] = y;
    } else {
      NumericVector z(l0.size());
      std::copy(y.begin(), y.end(), z.begin());
      lambda[i] = z;
    }
  }
  
  return lambda;
}
