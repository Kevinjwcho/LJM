#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
List calc_mu_surv(List X, List Y, List Z, List inv_omega, List K, arma::mat D, 
                  arma::vec beta, arma::vec phi, arma::vec eta, int nK) {
  int n = X.size();
  List mu_surv(n);
  
  for (int i = 0; i < n; i++) {
    arma::mat Xi = X[i];
    arma::vec Yi = Y[i];
    arma::mat Zi = Z[i];
    arma::mat inv_omegai = inv_omega[i];
    arma::mat Ki = K[i];
    
    // Calculate omega = diag(1/diag(inv_omega))
    arma::mat omega = inv_omegai;
    for (int j = 0; j < omega.n_rows; j++) {
      omega(j, j) = 1.0 / omega(j, j);
    }
    
    // Calculate mu_b = D * Z.t() * inv(omega + Z * D * Z.t()) * (Y - X * beta)
    arma::mat ZD = Zi * D; // Z * D
    arma::mat ZDZt = ZD * Zi.t();   // Z * D * Z.t()
    arma::mat term = omega + ZDZt;
    arma::mat term_inv = inv(term); // inv(omega + Z * D * Z.t())
    arma::vec Xbeta = Xi * beta; // X * beta
    arma::vec Y_Xbeta = Yi - Xbeta;
    arma::vec mu_b = ZD.t() * (term_inv * Y_Xbeta);
    
    // mu_surv = K * phi + eta.t() * mu_b
    arma::vec Kphi = Ki * phi;  // K * phi
    double etamu_b = arma::dot(eta, mu_b); // eta' * mu_b (dot product for vector multiplication)
    mu_surv[i] = Kphi + etamu_b;
  }
  
  return mu_surv;
}

// [[Rcpp::export]]
List calc_Sigma2_surv(List X, List Z, List inv_omega, arma::mat D, arma::vec eta) {
  int n = X.size();
  List Sigma2_surv(n);
  
  for (int i = 0; i < n; i++) {
    arma::mat Zi = Z[i];
    arma::mat inv_omegai = inv_omega[i];
    
    // Calculate omega = diag(1/diag(inv_omega))
    arma::mat omega = inv_omegai;
    for (int j = 0; j < omega.n_rows; j++) {
      omega(j, j) = 1.0 / omega(j, j);
    }
    
    // Sigma_b  = D - D * Z.t() * inv(omega + Z * D * Z.t()) * Z * D
    arma::mat ZD = Zi * D; // Z * D
    arma::mat ZDZt = ZD * Zi.t();   // Z * D * Z.t()
    arma::mat term = omega + ZDZt;
    arma::mat term_inv = inv(term); // inv(omega + Z * D * Z.t())
    arma::mat Sigma_b = D - (ZD.t() * (term_inv * (ZD)));
    
    // Sigma_surv = eta.t() * Sigma_b * eta
    double Sigma_surv_i = arma::as_scalar(eta.t() * Sigma_b * eta); // eta' * Sigma_b * eta
    Sigma2_surv[i] = Sigma_surv_i;
  }
  
  return Sigma2_surv;
}
