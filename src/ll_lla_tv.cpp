// ============================================================
// LLAJEL Time-Varying: Joint log-likelihood + gradient + Hessian (w.r.t b)
// Kernel-weighted longitudinal part + centered prior (b - c)
// + Time-varying association via etaBB matrix
//
// Exported to R:
//   ll_lla_tv(b, Y, Z, V, D, mi, K, Delta, l0i, l0u, phi, etaBB, w, c)
//   gradll_lla_tv(b, Y, Z, V, D, mi, K, Delta, l0i, l0u, phi, etaBB, w, c)
//   sdll_lla_tv(b, Z, D, V, K, l0u, phi, etaBB, w)
//
// Key difference from ll_lla:
//   - eta (q x 1 vector) is replaced by etaBB (nT x q matrix)
//     where nT = number of risk-set time points for this subject
//   - survival hazard: sum_t l0u[t] * exp(K*phi + etaBB[t,] * b)
//   - event term: Delta * (K*phi + etaBB[last,] * b)
// ============================================================

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// ------------------------------------------------------------
// Stable logdet for SPD matrices
// ------------------------------------------------------------
inline double logdet_spd_tv(const arma::mat& A) {
  double sign = 0.0, ld = 0.0;
  arma::log_det(ld, sign, A);
  if (sign <= 0) Rcpp::stop("Matrix is not positive definite (log_det sign <= 0).");
  return ld;
}

// ------------------------------------------------------------
// Robust SPD inverse + log-det.
// For univariate LLA fits the random-effects covariance D can be near-singular
// (intercept and slope columns collinear), so a plain inv_sympd / log_det
// throws.  We symmetrize, try inv_sympd, then ridge + retry, then pinv.
// All three D-using ops in this file (ll_lla_tv, gradll_lla_tv, sdll_lla_tv)
// route through this helper to share the same regularized D.
// ------------------------------------------------------------
inline void inv_logdet_spd_robust_tv(const arma::mat& D,
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
// ll_lla_tv: NEGATIVE joint log-likelihood for subject i (time-varying)
// [[Rcpp::export]]
double ll_lla_tv(
    const arma::vec& b,            // q x 1  random effects
    const arma::colvec& Y,         // mi x 1 longitudinal observations
    const arma::mat& Z,            // mi x q design matrix
    const arma::mat& V,            // mi x mi diagonal residual covariance
    const arma::mat& D,            // q x q  covariance matrix of random effects
    const int mi,                  // number of observations
    const arma::rowvec& K,         // 1 x ncw kernel covariate vector
    const int Delta,               // event indicator
    const double l0i,              // baseline hazard at event time
    const arma::rowvec& l0u,       // 1 x nT baseline hazard values at risk times
    const arma::vec& phi,          // ncw x 1 coefficients for K
    const arma::mat& etaBB,        // nT x q  time-varying association matrix (etaBBs[[i]])
    const arma::vec& w,            // mi x 1 kernel weights
    const arma::vec& c             // q x 1  LLA center
){
  // ---- input dimension checks ----
  if ((int)Y.n_elem != mi) Rcpp::stop("ll_lla_tv: length(Y) must equal mi.");
  if ((int)w.n_elem != mi) Rcpp::stop("ll_lla_tv: length(w) must equal mi.");
  if ((int)Z.n_rows != mi) Rcpp::stop("ll_lla_tv: nrow(Z) must equal mi.");
  if ((int)V.n_rows != mi || (int)V.n_cols != mi) Rcpp::stop("ll_lla_tv: V must be mi x mi.");
  if ((int)c.n_elem != (int)b.n_elem) Rcpp::stop("ll_lla_tv: length(c) must equal length(b).");
  if ((int)etaBB.n_cols != (int)b.n_elem) Rcpp::stop("ll_lla_tv: ncol(etaBB) must equal length(b).");

  const double eps = 1e-12;
  arma::vec wpos = arma::clamp(w, eps, datum::inf);

  // -------------------------------------------------------
  // (1) Longitudinal part: kernel-weighted (identical to ll_lla)
  // -------------------------------------------------------
  arma::colvec resid = Y - Z * b;
  arma::vec log_diag_V = arma::log(V.diag());
  double sum_w = arma::sum(wpos);
  double weighted_logdet = arma::dot(wpos, log_diag_V);
  arma::colvec Vinv_resid = arma::solve(V, resid, arma::solve_opts::fast);
  double quad = arma::as_scalar( (resid % wpos).t() * Vinv_resid );

  double ll_long =
    - 0.5 * sum_w * std::log(2.0 * M_PI)
    - 0.5 * weighted_logdet
    - 0.5 * quad;

  // -------------------------------------------------------
  // (2) Prior part: b ~ N(c, D) (identical to ll_lla)
  // -------------------------------------------------------
  int q = (int)b.n_elem;
  arma::vec bc = b - c;
  arma::mat D_inv;
  double logdetD;
  inv_logdet_spd_robust_tv(D, D_inv, logdetD);
  arma::colvec Dinv_bc = D_inv * bc;
  double quad_b = arma::as_scalar(bc.t() * Dinv_bc);

  double ll_b =
    - ( (double)q / 2.0 ) * std::log(2.0 * M_PI)
    - 0.5 * logdetD
    - 0.5 * quad_b;

  // -------------------------------------------------------
  // (3) Survival part: Time-varying Cox contribution
  //     etaBB is nT x q, each row = eta(t_j) for risk time j
  //     event term uses last row of etaBB (at event/censor time)
  //     hazard integral: sum_t l0u[t] * exp(K*phi + etaBB[t,]*b)
  // -------------------------------------------------------
  double Kphi = arma::as_scalar(K * phi);
  int nT = (int)l0u.n_elem;

  // Event term: Delta * (K*phi + etaBB[last,]*b)
  arma::rowvec etaBB_last = etaBB.row(nT - 1);   // last risk time
  double etab_event = arma::as_scalar(etaBB_last * b);
  double temp = (Delta == 1) ? std::log(l0i) : 0.0;

  // Cumulative hazard: sum_t l0u[t] * exp(K*phi + etaBB[t,]*b)
  double cum_haz = 0.0;
  for (int t = 0; t < nT; t++) {
    double etab_t = arma::as_scalar(etaBB.row(t) * b);
    cum_haz += l0u[t] * std::exp(Kphi + etab_t);
  }

  double ll_surv =
    temp
    + (double)Delta * (Kphi + etab_event)
    - cum_haz;

  return -1.0 * (ll_long + ll_b + ll_surv);
}


// ------------------------------------------------------------
// gradll_lla_tv: gradient of NEGATIVE log-likelihood w.r.t b (time-varying)
// [[Rcpp::export]]
arma::colvec gradll_lla_tv(
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
    const arma::mat& etaBB,
    const arma::vec& w,
    const arma::vec& c
){
  if ((int)Y.n_elem != mi) Rcpp::stop("gradll_lla_tv: length(Y) must equal mi.");
  if ((int)w.n_elem != mi) Rcpp::stop("gradll_lla_tv: length(w) must equal mi.");
  if ((int)Z.n_rows != mi) Rcpp::stop("gradll_lla_tv: nrow(Z) must equal mi.");
  if ((int)V.n_rows != mi || (int)V.n_cols != mi) Rcpp::stop("gradll_lla_tv: V must be mi x mi.");
  if ((int)c.n_elem != (int)b.n_elem) Rcpp::stop("gradll_lla_tv: length(c) must equal length(b).");
  if ((int)etaBB.n_cols != (int)b.n_elem) Rcpp::stop("gradll_lla_tv: ncol(etaBB) must equal length(b).");

  const double eps = 1e-12;
  arma::vec wpos = arma::clamp(w, eps, datum::inf);

  // Longitudinal gradient (same as ll_lla)
  arma::colvec resid = Y - Z * b;
  arma::colvec Vinv_resid = arma::solve(V, resid, arma::solve_opts::fast);
  arma::colvec A_resid = Vinv_resid % wpos;

  // Prior gradient (same as ll_lla); use robust D^{-1} (see helper at top)
  arma::mat D_inv_g;
  double logdetD_g;
  inv_logdet_spd_robust_tv(D, D_inv_g, logdetD_g);
  arma::colvec Dinv_bc = D_inv_g * (b - c);

  // Survival gradient (time-varying)
  double Kphi = arma::as_scalar(K * phi);
  int nT = (int)l0u.n_elem;

  // Event term gradient: Delta * etaBB[last,]'
  arma::rowvec etaBB_last = etaBB.row(nT - 1);
  arma::vec grad_event = (double)Delta * etaBB_last.t();

  // Cumulative hazard gradient: sum_t l0u[t] * exp(K*phi + etaBB[t,]*b) * etaBB[t,]'
  arma::vec grad_cumhaz = arma::zeros<arma::vec>(b.n_elem);
  for (int t = 0; t < nT; t++) {
    double etab_t = arma::as_scalar(etaBB.row(t) * b);
    double w_t = l0u[t] * std::exp(Kphi + etab_t);
    grad_cumhaz += w_t * etaBB.row(t).t();
  }

  // grad of LOG-likelihood wrt b:
  arma::colvec grad_loglik =
    Z.t() * A_resid     // longitudinal
    - Dinv_bc            // prior
    + grad_event         // survival event
    - grad_cumhaz;       // survival cumulative hazard

  return -1.0 * grad_loglik;
}


// ------------------------------------------------------------
// sdll_lla_tv: Hessian (2nd derivative) of LOG-likelihood wrt b (time-varying)
// [[Rcpp::export]]
arma::mat sdll_lla_tv(
    const arma::vec& b,
    const arma::mat& Z,
    const arma::mat& D,
    const arma::mat& V,
    const arma::rowvec& K,
    const arma::vec& l0u,
    const arma::vec& phi,
    const arma::mat& etaBB,
    const arma::vec& w
){
  int mi = (int)Z.n_rows;
  if ((int)w.n_elem != mi) Rcpp::stop("sdll_lla_tv: length(w) must equal nrow(Z).");
  if ((int)V.n_rows != mi || (int)V.n_cols != mi) Rcpp::stop("sdll_lla_tv: V must be mi x mi.");

  const double eps = 1e-12;
  arma::vec wpos = arma::clamp(w, eps, datum::inf);

  // Longitudinal Hessian: - Z^T diag(w) V^{-1} Z
  arma::mat VinvZ = arma::solve(V, Z, arma::solve_opts::fast);
  arma::mat A_Z = VinvZ.each_col() % wpos;
  arma::mat H_long = - (Z.t() * A_Z);

  // Prior Hessian: -D^{-1} (robust against near-singular D -- see helper)
  arma::mat D_inv;
  double logdetD_unused;
  inv_logdet_spd_robust_tv(D, D_inv, logdetD_unused);
  arma::mat H_prior = - D_inv;

  // Survival Hessian (time-varying):
  // - sum_t l0u[t] * exp(K*phi + etaBB[t,]*b) * etaBB[t,]' * etaBB[t,]
  double Kphi = arma::as_scalar(K * phi);
  int nT = (int)l0u.n_elem;
  int q = (int)b.n_elem;
  arma::mat H_surv = arma::zeros<arma::mat>(q, q);

  for (int t = 0; t < nT; t++) {
    double etab_t = arma::as_scalar(etaBB.row(t) * b);
    double kernel = l0u[t] * std::exp(Kphi + etab_t);
    H_surv -= kernel * (etaBB.row(t).t() * etaBB.row(t));
  }

  return H_long + H_prior + H_surv;
}
