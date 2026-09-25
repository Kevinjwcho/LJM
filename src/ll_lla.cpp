// ============================================================
// LLAJEL: Joint log-likelihood + gradient + Hessian (w.r.t b)
// Kernel-weighted longitudinal part + centered prior (b - c)
//
// Exported to R (names as requested):
//   ll_lla(b, Y, Z, V, D, mi, K, Delta, l0i, l0u, phi, eta, w, c)
//   gradll_lla(b, Y, Z, V, D, mi, K, Delta, l0i, l0u, phi, eta, w, c)
//   sdll_lla(b, Z, D, V, K, l0u, phi, eta, w)
//
// Conventions:
//   - ll_lla returns NEGATIVE log-likelihood (objective to minimize)
//   - gradll_lla returns gradient of the objective
//   - sdll_lla returns Hessian of the LOG-likelihood (same sign as old sdll).
//     If you want Hessian of the objective, multiply by (-1).
//
// Notes:
//   - V is subject-specific residual covariance (mi x mi), typically block-diag.
//   - w is kernel weights (length mi), aligned with stacked Y/Z rows.
//   - c is center vector from cLLA (length q), used in prior (b - c).
//   - Prefer filtering w>0 in R; we clamp to eps for safety.
//   - rho >= 0 is the transformation parameter of G(x) = log(1 + rho x)/rho
//     (rho = 0: Cox).  With u = Lambda0(V) exp(K phi + eta'b) the survival
//     log-density is Delta{log l0 + lp - log(1 + rho u)} - G(u).  rho = 0 takes
//     the original Cox code path unchanged, so results are bit-identical.
// ============================================================

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// ------------------------------------------------------------
// Stable logdet for SPD matrices
// ------------------------------------------------------------
inline double logdet_spd(const arma::mat& A) {
  double sign = 0.0, ld = 0.0;
  arma::log_det(ld, sign, A);
  if (sign <= 0) Rcpp::stop("Matrix is not positive definite (log_det sign <= 0).");
  return ld;
}

// ------------------------------------------------------------
// Robust SPD inverse + log-det (mirrors ll_lla_tv.cpp version).
// ------------------------------------------------------------
inline void inv_logdet_spd_robust(const arma::mat& D,
                                  arma::mat& D_inv,
                                  double& logdetD) {
  arma::mat D_sym = 0.5 * (D + D.t());
  bool ok = arma::inv_sympd(D_inv, D_sym);
  if (!ok) {
    double ridge = 1e-6 * arma::trace(D_sym)
      / static_cast<double>(std::max<arma::uword>(D_sym.n_rows, 1));
    if (ridge <= 0.0 || !std::isfinite(ridge)) ridge = 1e-8;
    D_sym = D_sym + ridge * arma::eye(D_sym.n_rows, D_sym.n_cols);
    ok = arma::inv_sympd(D_inv, D_sym);
    if (!ok) D_inv = arma::pinv(D_sym);
  }
  double sign = 0.0, ld = 0.0;
  arma::log_det(ld, sign, D_sym);
  logdetD = (sign > 0) ? ld : std::log(std::numeric_limits<double>::min());
}

// ------------------------------------------------------------
// ll_lla: NEGATIVE joint log-likelihood for subject i
// [[Rcpp::export]]
double ll_lla(
    const arma::vec& b,            // q x 1  random effects at landmark s
    const arma::colvec& Y,         // mi x 1 longitudinal observations
    const arma::mat& Z,            // mi x q design matrix: [1, t_ij - s]
    const arma::mat& V,            // mi x mi diagonal residual covariance diag(sigma_j^2)
    const arma::mat& D,            // q x q  covariance matrix of random effects
    const int mi,                  // number of observations for subject i
    const arma::rowvec& K,         // 1 x ncw kernel covariate vector (survival model)
    const int Delta,               // event indicator (1 = event, 0 = censored)
    const double l0i,              // baseline hazard at event time t_i
    const arma::rowvec& l0u,       // 1 x (#risk times) baseline hazard values at risk times
    const arma::vec& phi,          // ncw x 1 coefficients for K in survival model
    const arma::vec& eta,          // q x 1  coefficients for b in survival model
    const arma::vec& w,            // mi x 1 kernel weights: K_h(t_ij - s)
    const arma::vec& c,            // q x 1  LLA center (population-level intercept & slope at s)
    const double rho = 0.0         // transformation parameter (0 = Cox)
){
  // ---- input dimension checks ----
  if ((int)Y.n_elem != mi) Rcpp::stop("ll_lla: length(Y) must equal mi.");
  if ((int)w.n_elem != mi) Rcpp::stop("ll_lla: length(w) must equal mi.");
  if ((int)Z.n_rows != mi) Rcpp::stop("ll_lla: nrow(Z) must equal mi.");
  if ((int)V.n_rows != mi || (int)V.n_cols != mi) Rcpp::stop("ll_lla: V must be mi x mi.");
  if ((int)c.n_elem != (int)b.n_elem) Rcpp::stop("ll_lla: length(c) must equal length(b).");
  if ((int)eta.n_elem != (int)b.n_elem) Rcpp::stop("ll_lla: length(eta) must equal length(b).");
  
  // Clamp kernel weights to avoid log(0); prefer filtering w > 0 in R before calling
  const double eps = 1e-12;
  arma::vec wpos = arma::clamp(w, eps, datum::inf);
  
  // -------------------------------------------------------
  // (1) Longitudinal part: kernel-weighted log-likelihood
  //     ell_long = - sum_w/2 * log(2pi)
  //                - 1/2 * w' log(diag(V))
  //                - 1/2 * (Y-Zb)' W V^{-1} (Y-Zb)
  //
  //     V is diagonal: V = diag(sigma_1^2, ..., sigma_m^2)
  //     W = diag(w_1, ..., w_m), w_j = K_h(t_ij - s)
  // -------------------------------------------------------
  
  // Residual: r = Y - Z*b
  arma::colvec resid = Y - Z * b;
  
  // log(sigma_j^2) = log of diagonal elements of V
  arma::vec log_diag_V = arma::log(V.diag());
  
  // sum of kernel weights: sum_j w_j
  double sum_w = arma::sum(wpos);
  
  // kernel-weighted log-det term: sum_j w_j * log(sigma_j^2)
  double weighted_logdet = arma::dot(wpos, log_diag_V);
  
  // quadratic form: r' W V^{-1} r
  // V is diagonal so V^{-1} r = r / diag(V)
  arma::colvec Vinv_resid = arma::solve(V, resid, arma::solve_opts::fast);
  double quad = arma::as_scalar( (resid % wpos).t() * Vinv_resid );
  
  double ll_long =
    - 0.5 * sum_w * std::log(2.0 * M_PI)  // - sum(w)/2 * log(2pi)
    - 0.5 * weighted_logdet                // - 1/2 * sum_j w_j * log(sigma_j^2)
    - 0.5 * quad;                          // - 1/2 * r' W V^{-1} r
    
    // -------------------------------------------------------
    // (2) Prior part: b ~ N(c, D)
    //     c = LLA population-level center (intercept & slope at s)
    //     ell_prior = -q/2 * log(2pi) - 1/2 * log|D|
    //                 - 1/2 * (b-c)' D^{-1} (b-c)
    // -------------------------------------------------------
    int q = (int)b.n_elem;
    arma::vec bc = b - c;  // centered random effect

    arma::mat D_inv;
    double logdetD;
    inv_logdet_spd_robust(D, D_inv, logdetD);
    arma::colvec Dinv_bc = D_inv * bc;
    double quad_b = arma::as_scalar(bc.t() * Dinv_bc);
    
    double ll_b =
      - ( (double)q / 2.0 ) * std::log(2.0 * M_PI)
      - 0.5 * logdetD
      - 0.5 * quad_b;
      
      // -------------------------------------------------------
      // (3) Survival part: Cox proportional hazards contribution
      //     linear predictor = K*phi + eta'*b
      //       K*phi  : kernel covariate term (time-fixed at landmark s)
      //       eta'*b : random effect term
      //     ell_surv = Delta * log(l0i)
      //                + Delta * (K*phi + eta'*b)
      //                - sum(l0u) * exp(K*phi) * exp(eta'*b)
      // -------------------------------------------------------
      double Kphi  = arma::as_scalar(K * phi);    // scalar: K * phi
      double etab  = arma::as_scalar(eta.t() * b); // scalar: eta' * b
      
      double temp     = (Delta == 1) ? std::log(l0i) : 0.0;  // log(l0i) only if event occurred
      double sum_l0u  = arma::sum(l0u);                        // cumulative baseline hazard
      
      double ll_surv;
      if (rho == 0.0) {
        ll_surv =
        temp
        + (double)Delta * (Kphi + etab)                    // Delta * linear predictor
        - sum_l0u * std::exp(Kphi) * std::exp(etab);       // expected cumulative hazard
      } else {
        // Transformation model: Delta{log l0 + lp + log G'(u)} - G(u),
        // G'(u) = 1/(1 + rho u),  G(u) = log(1 + rho u)/rho.
        double u      = sum_l0u * std::exp(Kphi) * std::exp(etab);
        double log1ru = std::log1p(rho * u);
        ll_surv = temp
          + (double)Delta * (Kphi + etab)
          - (double)Delta * log1ru
          - log1ru / rho;
      }
        
        // -------------------------------------------------------
        // Return NEGATIVE joint pseudo-log-likelihood (to minimize)
        // -------------------------------------------------------
        return -1.0 * (ll_long + ll_b + ll_surv);
}


// ------------------------------------------------------------
// gradll_lla: gradient of NEGATIVE log-likelihood w.r.t b
// [[Rcpp::export]]
arma::colvec gradll_lla(
    const arma::vec& b,
    const arma::colvec& Y,
    const arma::mat& Z,
    const arma::mat& V,
    const arma::mat& D,
    const int mi,
    const arma::rowvec& K,
    const int Delta,
    const double l0i,
    const arma::rowvec& l0u,
    const arma::vec& phi,
    const arma::vec& eta,
    const arma::vec& w,
    const arma::vec& c,
    const double rho = 0.0
){
  // ---- minimal checks ----
  if ((int)Y.n_elem != mi) Rcpp::stop("gradll_lla: length(Y) must equal mi.");
  if ((int)w.n_elem != mi) Rcpp::stop("gradll_lla: length(w) must equal mi.");
  if ((int)Z.n_rows != mi) Rcpp::stop("gradll_lla: nrow(Z) must equal mi.");
  if ((int)V.n_rows != mi || (int)V.n_cols != mi) Rcpp::stop("gradll_lla: V must be mi x mi.");
  if ((int)c.n_elem != (int)b.n_elem) Rcpp::stop("gradll_lla: length(c) must equal length(b).");
  if ((int)eta.n_elem != (int)b.n_elem) Rcpp::stop("gradll_lla: length(eta) must equal length(b).");
  
  const double eps = 1e-12;
  arma::vec wpos = arma::clamp(w, eps, datum::inf);
  
  // residual
  arma::colvec resid = Y - Z * b;
  
  // A_resid = diag(w) V^{-1} resid
  arma::colvec Vinv_resid = arma::solve(V, resid, arma::solve_opts::fast);
  arma::colvec A_resid    = Vinv_resid % wpos;
  
  // survival hazard multiplier
  double sum_l0u = arma::sum(l0u);
  double expKphi = std::exp( arma::as_scalar(K * phi) );
  double expetab = std::exp( arma::as_scalar(eta.t() * b) );
  double haz_mult = sum_l0u * expKphi * expetab;
  
  arma::vec Delta_vec = (double)Delta * arma::ones<vec>(eta.n_elem);
  
  // grad of LOG-likelihood wrt b (use robust D^{-1} — see helper at top):
  arma::mat D_inv_g;
  double logdetD_g;
  inv_logdet_spd_robust(D, D_inv_g, logdetD_g);

  arma::colvec grad_loglik;
  if (rho == 0.0) {
    grad_loglik =
    Z.t() * A_resid
  - D_inv_g * (b - c)
    + Delta_vec % eta
    - haz_mult * eta;
  } else {
    // d/db [Delta{lp - log(1+rho u)} - G(u)] = eta (Delta - u)/(1 + rho u)
    double gsurv = ((double)Delta - haz_mult) / (1.0 + rho * haz_mult);
    grad_loglik =
      Z.t() * A_resid
      - D_inv_g * (b - c)
      + gsurv * eta;
  }
    
    // gradient of NEGATIVE log-likelihood
    return -1.0 * grad_loglik;
}


// ------------------------------------------------------------
// sdll_lla: Hessian (2nd derivative) of LOG-likelihood wrt b
// (same sign convention as your old sdll)
// [[Rcpp::export]]
arma::mat sdll_lla(
    const arma::vec& b,
    const arma::mat& Z,
    const arma::mat& D,
    const arma::mat& V,
    const arma::rowvec& K,
    const arma::vec& l0u,
    const arma::vec& phi,
    const arma::vec& eta,
    const arma::vec& w,
    const int Delta = 0,           // needed only when rho > 0
    const double rho = 0.0
){
  int mi = (int)Z.n_rows;
  if ((int)w.n_elem != mi) Rcpp::stop("sdll_lla: length(w) must equal nrow(Z).");
  if ((int)V.n_rows != mi || (int)V.n_cols != mi) Rcpp::stop("sdll_lla: V must be mi x mi.");
  
  const double eps = 1e-12;
  arma::vec wpos = arma::clamp(w, eps, datum::inf);
  
  // Longitudinal Hessian:  - Z^T diag(w) V^{-1} Z
  arma::mat VinvZ = arma::solve(V, Z, arma::solve_opts::fast); // mi x q
  arma::mat A_Z   = VinvZ.each_col() % wpos;                  // diag(w) * (V^{-1} Z)
  arma::mat H_long = - (Z.t() * A_Z);
  
  // Prior Hessian: -D^{-1} (robust against near-singular D — see helper)
  arma::mat D_inv;
  double logdetD_unused;
  inv_logdet_spd_robust(D, D_inv, logdetD_unused);
  arma::mat H_prior = - D_inv;
  
  // Survival Hessian: -kernel*(eta eta^T)
  double sum_l0u = arma::sum(l0u);
  double kernel  = sum_l0u
  * std::exp( arma::as_scalar(K * phi) )
    * std::exp( arma::as_scalar(eta.t() * b) );
    
    arma::mat H_surv;
    if (rho == 0.0) {
      H_surv = - kernel * (eta * eta.t());
    } else {
      // d2/db2 = - u (1 + rho Delta) / (1 + rho u)^2  eta eta'
      double opr = 1.0 + rho * kernel;
      H_surv = - (kernel * (1.0 + rho * (double)Delta) / (opr * opr)) * (eta * eta.t());
    }
    
    return H_long + H_prior + H_surv;
}