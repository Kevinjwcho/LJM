#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List update_lambda(NumericVector Es_exp, List l0u, NumericVector l0, int n) {
  List lambda(n);
  
  for (int i = 0; i < n; i++) {
    NumericVector l0u_i = l0u[i];
    int len_l0u = l0u_i.size();
    
    NumericVector y(len_l0u);
    std::fill(y.begin(), y.end(), Es_exp[i]);
    
    if (len_l0u == l0.size()) {
      lambda[i] = y;
    } else {
      NumericVector z(l0.size());
      std::copy(y.begin(), y.end(), z.begin());
      lambda[i] = z;
    }
  }
  
  return lambda;
}
