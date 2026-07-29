#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Joint log-likelihood multiplying minus to minimize
// [[Rcpp::export]]
double ll(arma::vec b, const arma::colvec Y, const arma::mat X, const arma::mat Z, const arma::mat V, const arma::mat D,
		   int mi, const arma::rowvec K, const int Delta, const double l0i,
		   const arma::rowvec l0u,
		   const arma::vec beta, const arma::vec phi,
		   const arma::vec eta,
		   const int nK, const int q){
	double temp = 0.0;
	if(Delta == 1) temp = log(l0i);
	// Sum of l0u
	double sum_l0u = arma::sum(l0u);
	arma::colvec resid = Y - X * beta - Z * b;
	return -1.0 * as_scalar(-mi/2.0 * log(2.0 * M_PI) - 0.5 * log(det(V)) -0.5 * resid.t() * V.i() * resid +
	                 -nK * log(2.0 * M_PI) - 0.5 * log(det(D)) - 0.5 * b.t() * D.i() * b +
	                 temp + Delta * (K * phi + (eta.t() * b)) - sum_l0u * (exp(K * phi) * exp((eta.t() * b))));
}

// Gradient function, NB some redundant arguments included for ease of use through ucminf with ll() above multipying minus
// [[Rcpp::export]]
arma::colvec gradll(arma::vec b, const arma::colvec Y, const arma::mat X, const arma::mat Z, const arma::mat V, const arma::mat D,
		     int mi, const arma::rowvec K, const int Delta, const double l0i,
		     const arma::rowvec l0u,
		     const arma::vec beta, const arma::vec phi,
		     const arma::vec eta,
		     const int nK, const int q){
  arma::colvec resid = Y - X * beta - Z * b;
  arma::vec Delta_vec = Delta * ones<vec>(eta.n_elem);  // Creating a vector with the same size as eta filled with Delta
  double sum_l0u = arma::sum(l0u);
// 	return -1.0 * (Z.t() * V.i() * resid - D.i() * b + Delta_vec % eta +
// 	                 - ((sum_l0u * exp(K * phi) * exp((eta.t() * b))) * ones<vec>(eta.n_elem)) % eta
// );
  arma::colvec exp_term = exp(as_scalar(K * phi)) * exp(as_scalar(eta.t() * b)) * ones<vec>(eta.n_elem);  // Ensure scalar multiplication
  return -1.0 * (Z.t() * V.i() * resid - D.i() * b + Delta_vec % eta + (-sum_l0u * exp_term) % eta);
  return exp_term;

}

// Joint likelihood: second derivative w.r.t random effects b_i
// [[Rcpp::export]]
arma::mat sdll(arma::vec b, const arma::mat Z, const arma::mat D, const arma::mat V,
		 const arma::rowvec K, const arma::vec l0u, const arma::vec phi,
		 const arma::vec eta, const int nK){
  double sum_l0u = arma::sum(l0u);
  double kernel = as_scalar(sum_l0u * exp(K * phi) * exp(eta.t()*b));
  // mat diagkern = diagmat(kernel);
  return -1.0 * (Z.t() * V.i() * Z) - D.i()
    + (-kernel*(eta * eta.t()))
    ;
}

