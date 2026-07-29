#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//////////////////////////////////////////////////////////////
// mu_surv (NO X, NO beta)
// \tilde Y = Y - Zc
// mu_b = c + D * Z^T * (omega + Z D Z^T)^{-1} * \tilde Y
// mu_surv = K * phi + (eta^T * mu_b) * 1
//
// Assumptions:
// - inv_omega is diagonal precision matrix (kernel-weighted): diag(w / var_e)
// - omega is diagonal variance matrix: diag(var_e / w) = diag(1 / diag(inv_omega))
// - K_i * phi returns a vector (same length as what you want for mu_surv_i)
// - eta is a vector compatible with mu_b (same length as random effects dimension)
//////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List calc_mu_surv_noX(List Y, List Z, List inv_omega, List K,
                      const arma::mat& D,
                      const arma::vec& phi,
                      const arma::vec& eta,
                      List c){  // q x 1, LLA center: c = (c0, c1)
  // c0: population-level intercept at s
  // c1: population-level slope at s
  int n = Y.size();
  List mu_surv(n);
  
  for (int i = 0; i < n; i++) {
    arma::vec Yi         = as<arma::vec>(Y[i]);
    arma::mat Zi         = as<arma::mat>(Z[i]);  // mi x q: [1, t_ij - s]
    arma::mat inv_omegai = as<arma::mat>(inv_omega[i]);
    arma::mat Ki         = as<arma::mat>(K[i]);
    arma::vec ci         = as<arma::vec>(c[i]); 
    
    // omega = diag(1 / diag(inv_omega))
    arma::vec prec = inv_omegai.diag();
    if (prec.min() <= 0) stop("inv_omega diagonal must be positive (subject %d).", i + 1);
    arma::vec omega_diag = 1.0 / prec;
    
    // -------------------------------------------------------
    // Center Y by population-level trajectory:
    //   Z_i * c = [c0 + c1*(t_i1-s), ..., c0 + c1*(t_im-s)]
    //           = population-level trajectory at observed times
    //
    // Y_tilde = Y_i - Z_i * c  (individual deviation from population)
    // This shifts b ~ N(c, D) to b_tilde ~ N(0, D)
    // -------------------------------------------------------
    arma::vec Zic      = Zi * ci;           // mi x 1: population trajectory
    arma::vec Yi_tilde = Yi - Zic;         // mi x 1: individual deviation
    
    // term = Omega_i + Z_i D Z_i^T
    arma::mat ZD   = Zi * D;
    arma::mat term = arma::diagmat(omega_diag) + (ZD * Zi.t());
    
    // Solve term * x = Y_tilde
    arma::vec x = arma::solve(term, Yi_tilde, arma::solve_opts::fast);
    
    // BLUP of b_tilde = D Z_i^T x  (centered, b_tilde ~ N(0, D))
    arma::vec b_tilde = ZD.t() * x;
    
    // -------------------------------------------------------
    // Restore original scale: b_hat = b_tilde + c
    //   b_tilde: individual deviation from population (intercept & slope)
    //   c      : population-level intercept & slope at landmark s
    //   b_hat  : BLUP of b_i ~ N(c, D)
    // -------------------------------------------------------
    arma::vec b_hat = b_tilde + ci;        // q x 1: (b0_hat, b1_hat)
    
    // mu_surv_i = K_i * phi + eta' * b_hat
    //   K_i * phi : kernel covariate term (time-fixed at landmark s)
    //   eta'*b_hat: random effect contribution to survival
    arma::vec Kphi   = Ki * phi;
    double eta_b_hat = arma::dot(eta, b_hat);
    
    mu_surv[i] = Kphi + eta_b_hat;
  }
  
  return mu_surv;
}

//////////////////////////////////////////////////////////////
// Sigma2_surv (NO X)
// Sigma_b = D - (Z D)^T * (omega + Z D Z^T)^{-1} * (Z D)
// Sigma2_surv = eta^T Sigma_b eta  (scalar)
// Kernel enters through inv_omega -> omega
//////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List calc_Sigma2_surv_noX(List Z, List inv_omega,
                          const arma::mat& D,
                          const arma::vec& eta) {
  int n = Z.size();
  List Sigma2_surv(n);
  
  for (int i = 0; i < n; i++) {
    arma::mat Zi = as<arma::mat>(Z[i]);
    arma::mat inv_omegai = as<arma::mat>(inv_omega[i]);
    
    arma::vec prec = inv_omegai.diag();
    if (prec.min() <= 0) stop("inv_omega diagonal must be positive (subject %d).", i + 1);
    arma::vec omega_diag = 1.0 / prec;
    
    arma::mat ZD = Zi * D;
    arma::mat term = arma::diagmat(omega_diag) + (ZD * Zi.t());
    
    // Compute A = term^{-1} * (Z D) via solve
    arma::mat A = arma::solve(term, ZD, arma::solve_opts::fast);
    
    // Sigma_b = D - (Z D)^T * A
    arma::mat Sigma_b = D - (ZD.t() * A);
    
    double Sigma_surv_i = arma::as_scalar(eta.t() * Sigma_b * eta);
    Sigma2_surv[i] = Sigma_surv_i;
  }
  
  return Sigma2_surv;
}
