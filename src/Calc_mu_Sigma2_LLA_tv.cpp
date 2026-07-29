#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Forward declarations for functions defined in etaphiRcpp_tv.cpp
List etaBBs_compute(const List& BBi, const arma::vec& eta, const arma::vec& etaInd);
arma::vec Esurv_tv(arma::vec w, arma::vec v, List mu, List variance,
                   List mu_new, List variance_new, arma::vec l0i, List l0u);
arma::vec Esurv_exp_tv(arma::vec w, arma::vec v, List mu, List variance,
                       List mu_new, List variance_new, arma::vec l0i, List l0u);

//////////////////////////////////////////////////////////////
// Time-varying mu_surv (NO X, NO beta, with etaBB matrix)
//
// mu_surv_i = Kron(K*phi, 1_{nT}) + etaBB * mu_b
//   where mu_b = c + D * Z^T * (omega + Z D Z^T)^{-1} * (Y - Z*c)
//   and etaBB is nT x q (time-varying association)
//
// Returns a vector of length nT for each subject (not scalar)
//////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List calc_mu_surv_noX_tv(List Y, List Z, List inv_omega, List K,
                         const arma::mat& D,
                         const arma::vec& phi,
                         List etaBB,       // List of nT_i x q matrices
                         List c){
  int n = Y.size();
  List mu_surv(n);

  for (int i = 0; i < n; i++) {
    arma::vec Yi         = as<arma::vec>(Y[i]);
    arma::mat Zi         = as<arma::mat>(Z[i]);
    arma::mat inv_omegai = as<arma::mat>(inv_omega[i]);
    arma::mat Ki         = as<arma::mat>(K[i]);
    arma::vec ci         = as<arma::vec>(c[i]);
    arma::mat etaB       = as<arma::mat>(etaBB[i]);  // nT_i x q

    // omega = diag(1 / diag(inv_omega))
    arma::vec prec = inv_omegai.diag();
    if (prec.min() <= 0) stop("inv_omega diagonal must be positive (subject %d).", i + 1);
    arma::vec omega_diag = 1.0 / prec;

    // Center Y by population trajectory: Y_tilde = Y - Z*c
    arma::vec Zic      = Zi * ci;
    arma::vec Yi_tilde = Yi - Zic;

    // term = Omega + Z D Z^T
    arma::mat ZD   = Zi * D;
    arma::mat term = arma::diagmat(omega_diag) + (ZD * Zi.t());

    // BLUP of b_tilde = D Z^T (term)^{-1} Y_tilde
    arma::vec x = arma::solve(term, Yi_tilde, arma::solve_opts::fast);
    arma::vec b_tilde = ZD.t() * x;

    // Restore: b_hat = b_tilde + c
    arma::vec b_hat = b_tilde + ci;

    // mu_surv_i = kron(K*phi, ones(nT)) + etaBB * b_hat
    arma::vec Kphi = arma::kron((Ki * phi), arma::ones<arma::vec>(etaB.n_rows));
    arma::vec etamu_b = etaB * b_hat;   // nT x 1

    mu_surv[i] = Kphi + etamu_b;
  }

  return mu_surv;
}

//////////////////////////////////////////////////////////////
// Time-varying Sigma2_surv (NO X, with etaBB matrix)
//
// Sigma_b = D - (Z D)^T * (omega + Z D Z^T)^{-1} * (Z D)
// Sigma2_surv = diag(etaBB * Sigma_b * etaBB^T)   (vector of length nT)
//////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List calc_Sigma2_surv_noX_tv(List Z, List inv_omega,
                             const arma::mat& D,
                             List etaBB) {
  int n = Z.size();
  List Sigma2_surv(n);

  for (int i = 0; i < n; i++) {
    arma::mat Zi = as<arma::mat>(Z[i]);
    arma::mat inv_omegai = as<arma::mat>(inv_omega[i]);
    arma::mat etaB = as<arma::mat>(etaBB[i]);  // nT_i x q

    arma::vec prec = inv_omegai.diag();
    if (prec.min() <= 0) stop("inv_omega diagonal must be positive (subject %d).", i + 1);
    arma::vec omega_diag = 1.0 / prec;

    arma::mat ZD = Zi * D;
    arma::mat term = arma::diagmat(omega_diag) + (ZD * Zi.t());

    // A = term^{-1} * (Z D)
    arma::mat A = arma::solve(term, ZD, arma::solve_opts::fast);

    // Sigma_b = D - (Z D)^T * A
    arma::mat Sigma_b = D - (ZD.t() * A);

    // Sigma_surv_i = etaBB * Sigma_b * etaBB^T  ->  diag()
    arma::mat Sigma_surv_i = etaB * Sigma_b * etaB.t();
    arma::vec diag_sigma_surv = arma::vec(Sigma_surv_i.diag());

    Sigma2_surv[i] = wrap(diag_sigma_surv);
  }

  return Sigma2_surv;
}


//////////////////////////////////////////////////////////////
// Eetaphi_lla_tv: Expected log-likelihood for (eta, phi) update
// LLA version (no X, no beta, kernel-weighted posterior, centered prior)
//
// This parallels Eetaphi_t from etaphiRcpp_tv.cpp but:
//   - No X/beta in mu_b computation
//   - Uses kernel-weighted inv_omega (inv_omega_ker)
//   - Uses centered prior: mu_b = c + D Z^T (omega + Z D Z^T)^{-1} (Y - Zc)
//////////////////////////////////////////////////////////////

// [[Rcpp::export]]
double Eetaphi_lla_tv(arma::vec etaphi, List Y, List Z, List inv_omega, List K,
                      arma::mat D, arma::vec l0i, List l0u,
                      arma::vec Di, int nK, arma::vec w, arma::vec v,
                      List mu_surv, List Sigma2_surv,
                      const List& BBi, const arma::vec& etaInd,
                      List c_list) {

  int n_etaphi = etaphi.n_elem;
  int n_eta = sum(etaInd);
  arma::vec e = etaphi.subvec(0, n_eta - 1);
  arma::vec p = etaphi.subvec(n_eta, n_etaphi - 1);

  List eBBs = etaBBs_compute(BBi, e, etaInd);

  List mu_surv_new(Y.size());
  List Sigma2_surv_new(Y.size());

  for (int i = 0; i < Y.size(); i++) {
    arma::vec Yi         = as<arma::vec>(Y[i]);
    arma::mat Zi         = as<arma::mat>(Z[i]);
    arma::mat inv_omegai = as<arma::mat>(inv_omega[i]);
    arma::mat Ki         = as<arma::mat>(K[i]);
    arma::vec ci         = as<arma::vec>(c_list[i]);
    arma::mat eBBi       = as<arma::mat>(eBBs[i]);

    // omega = diag(1 / diag(inv_omega))
    arma::vec prec = inv_omegai.diag();
    arma::vec omega_diag = 1.0 / prec;

    // Centered BLUP: b_hat = c + D Z^T (Omega + Z D Z^T)^{-1} (Y - Zc)
    arma::vec Yi_tilde = Yi - Zi * ci;
    arma::mat ZD = Zi * D;
    arma::mat V_mat = arma::diagmat(omega_diag) + ZD * Zi.t();
    arma::vec mu_b = ci + ZD.t() * arma::solve(V_mat, Yi_tilde, arma::solve_opts::fast);

    // mu_surv_i = kron(K*p, ones) + eBBi * mu_b
    arma::vec mu_surv_i = arma::ones<arma::vec>(eBBi.n_rows) * arma::as_scalar(Ki * p) + eBBi * mu_b;
    mu_surv_new[i] = mu_surv_i;

    // Sigma_b = D - D Z^T (V)^{-1} Z D
    arma::mat Sigma_b = D - ZD.t() * arma::solve(V_mat, ZD, arma::solve_opts::fast);
    arma::mat Sigma_surv_i = eBBi * Sigma_b * eBBi.t();
    arma::vec diag_sigma = arma::vec(Sigma_surv_i.diag());
    Sigma2_surv_new[i] = wrap(diag_sigma);
  }

  // Re-use Esurv_tv and Esurv_exp_tv from etaphiRcpp_tv.cpp (already compiled)
  arma::vec Es_exp = Esurv_exp_tv(w, v, mu_surv, Sigma2_surv,
                                  mu_surv_new, Sigma2_surv_new, l0i, l0u);
  arma::vec Es = Esurv_tv(w, v, mu_surv, Sigma2_surv,
                          mu_surv_new, Sigma2_surv_new, l0i, l0u);

  double result = 0.0;
  for (unsigned int i = 0; i < Di.n_elem; i++) {
    result += Di[i] * Es[i] - Es_exp[i];
  }

  return result;
}


// [[Rcpp::export]]
arma::vec Setaphi_lla_tv(arma::vec etaphi,
                         List Y, List Z, List inv_omega,
                         List K, arma::mat D,
                         arma::vec l0i, List l0u, arma::vec Di,
                         int nK, arma::vec w, arma::vec v,
                         List mu_surv, List Sigma2_surv,
                         const List& BBi, const arma::vec& etaInd,
                         double eps, List c_list) {

  int ep_size = etaphi.n_elem;
  arma::vec gradient(ep_size);
  double f0 = Eetaphi_lla_tv(etaphi, Y, Z, inv_omega, K, D, l0i, l0u, Di, nK, w, v,
                             mu_surv, Sigma2_surv, BBi, etaInd, c_list);

  for (int i = 0; i < ep_size; ++i) {
    arma::vec ep = etaphi;
    double xi = std::max(std::abs(ep[i]), 1.0);
    double h = xi * eps;
    ep[i] += h;
    double f1 = Eetaphi_lla_tv(ep, Y, Z, inv_omega, K, D, l0i, l0u, Di, nK, w, v,
                               mu_surv, Sigma2_surv, BBi, etaInd, c_list);
    gradient[i] = (f1 - f0) / h;
  }

  return gradient;
}


// [[Rcpp::export]]
arma::mat Hetaphi_lla_tv(arma::vec etaphi,
                         List Y, List Z, List inv_omega,
                         List K, arma::mat D,
                         arma::vec l0i, List l0u, arma::vec Di,
                         int nK, arma::vec w, arma::vec v,
                         List mu_surv, List Sigma2_surv,
                         const List& BBi, const arma::vec& etaInd,
                         double eps, List c_list) {

  int ep_size = etaphi.n_elem;
  arma::mat hessian(ep_size, ep_size);
  arma::vec f0 = Setaphi_lla_tv(etaphi, Y, Z, inv_omega, K, D,
                                l0i, l0u, Di, nK, w, v,
                                mu_surv, Sigma2_surv, BBi, etaInd, eps, c_list);

  for (int i = 0; i < ep_size; ++i) {
    arma::vec ep = etaphi;
    double xi = std::max(std::abs(ep[i]), 1.0);
    double h = xi * eps;
    ep[i] += h;
    arma::vec f1 = Setaphi_lla_tv(ep, Y, Z, inv_omega, K, D,
                                  l0i, l0u, Di, nK, w, v,
                                  mu_surv, Sigma2_surv, BBi, etaInd, eps, c_list);
    hessian.col(i) = (f1 - f0) / h;
  }

  hessian = 0.5 * (hessian + hessian.t());

  return hessian;
}
