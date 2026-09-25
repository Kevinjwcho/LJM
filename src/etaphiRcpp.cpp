// Filename: etaphiRcpp.cpp  (NO X, NO beta; kernel-aware via inv_omega)
//
// ----------------------------------------------------------------------------
// Optimization note (2026-04):
// The original Eetaphi redundantly recomputed per-subject matrix work
// (V, ZD, solve(V,.), mu_b, Sigma_b) that does NOT depend on etaphi on every
// call. Inside Hetaphi, this redundancy compounded to (q+1)^2 full rebuilds,
// giving K^~2.87 scaling and >80% of EM wall time at K=5.
//
// Refactor: split Eetaphi into
//   (1) build_precomp(...)  : etaphi-independent per-subject (mu_b, Sigma_b)
//   (2) Eetaphi_fast(...)   : etaphi-dependent inner loop (3 scalars/subject)
// The exported Eetaphi/Setaphi/Hetaphi now build precomp ONCE per call and
// reuse it for all finite-difference evaluations. Numerical results are
// bit-identical to the previous implementation.
// ----------------------------------------------------------------------------

#include <RcppArmadillo.h>
#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//////////////////////////////////////////////////////////////
// Transformation model (rho >= 0, G(x) = log(1 + rho x)/rho; rho = 0: Cox).
// The 1-D quadrature weights use the survival density
//   f(V, Delta | zeta) = {l0 e^zeta G'(u)}^Delta exp{-G(u)},  u = Lambda0(V) e^zeta,
// and Esurv_exp returns E[xi e^{zeta_new} | O] with
//   E(xi | b, O) = (1 + Delta rho)/(1 + rho u)   (Zeng & Lin 2007).
// Delta is read from l0i (l0i == 0 <=> censored), as before.  rho = 0 takes
// the original Cox code path unchanged.
//////////////////////////////////////////////////////////////
// Esurv / Esurv_exp: 그대로 사용 (mu[[i]][0], var[[i]][0]만 사용)
//////////////////////////////////////////////////////////////

// [[Rcpp::export]]
arma::vec Esurv(arma::vec w, arma::vec v, List mu, List variance,
                List mu_new, List variance_new, arma::vec l0i, List l0u,
                double rho = 0.0) {

  int n = mu.size();
  int nq = w.n_elem;   // number of quadrature nodes
  arma::vec nom(n, fill::zeros);
  arma::vec denom(n, fill::zeros);

  for (int l = 0; l < nq; l++) {
    double ww = w[l];
    double vv = v[l];

    for (int i = 0; i < n; i++) {
      arma::vec mui = mu[i];
      arma::vec variancei = variance[i];
      arma::vec mu_newi = mu_new[i];
      arma::vec variance_newi = variance_new[i];

      double mui_val = mui[0];
      double variancei_val = variancei[0];
      double mu_newi_val = mu_newi[0];
      double variance_newi_val = variance_newi[0];

      double exp_mu = std::exp(mui_val + std::sqrt(variancei_val) * vv);
      double mu_new_val = mu_newi_val + std::sqrt(variance_newi_val) * vv;

      double firstTerm = (l0i[i] == 0) ? 1.0 : exp_mu * l0i[i];

      arma::vec l0u_i = l0u[i];
      double sum_l0u_exp = arma::sum(l0u_i * exp_mu);
      double secondTerm = std::exp(-sum_l0u_exp);
      if (rho != 0.0) {
        // G-model survival density: l0 e^zeta G'(u) for events, exp(-G(u))
        double log1ru = std::log1p(rho * sum_l0u_exp);
        if (l0i[i] != 0) firstTerm = exp_mu * l0i[i] / (1.0 + rho * sum_l0u_exp);
        secondTerm = std::exp(-log1ru / rho);
      }

      double denom_increment = ww * firstTerm * secondTerm;
      denom[i] += denom_increment;
      nom[i] += ww * mu_new_val * firstTerm * secondTerm;
    }
  }

  arma::vec result(n);
  for (int i = 0; i < n; i++) {
    result[i] = (denom[i] != 0) ? (nom[i] / denom[i]) : datum::nan;
  }
  return result;
}

// [[Rcpp::export]]
arma::vec Esurv_exp(arma::vec w, arma::vec v, List mu, List variance,
                    List mu_new, List variance_new, arma::vec l0i, List l0u,
                    double rho = 0.0) {

  int n = mu.size();
  int nq = w.n_elem;   // number of quadrature nodes
  arma::vec nom(n, fill::zeros);
  arma::vec denom(n, fill::zeros);

  for (int l = 0; l < nq; l++) {
    double ww = w[l];
    double vv = v[l];

    for (int i = 0; i < n; i++) {
      arma::vec mui = mu[i];
      arma::vec variancei = variance[i];
      arma::vec mu_newi = mu_new[i];
      arma::vec variance_newi = variance_new[i];

      double mui_val = mui[0];
      double variancei_val = variancei[0];
      double mu_newi_val = mu_newi[0];
      double variance_newi_val = variance_newi[0];

      double exp_mu = std::exp(mui_val + std::sqrt(variancei_val) * vv);
      double exp_mu_new = std::exp(mu_newi_val + std::sqrt(variance_newi_val) * vv);

      double firstTerm = (l0i[i] == 0) ? 1.0 : exp_mu * l0i[i];

      arma::vec l0u_i = l0u[i];
      double sum_l0u_exp = arma::sum(l0u_i * exp_mu);
      double secondTerm = std::exp(-sum_l0u_exp);

      if (rho == 0.0) {
        double denom_increment = ww * firstTerm * secondTerm;
        denom[i] += denom_increment;
        nom[i] += ww * exp_mu_new * firstTerm * secondTerm;
      } else {
        double log1ru = std::log1p(rho * sum_l0u_exp);
        double delta_i = (l0i[i] != 0) ? 1.0 : 0.0;
        if (l0i[i] != 0) firstTerm = exp_mu * l0i[i] / (1.0 + rho * sum_l0u_exp);
        secondTerm = std::exp(-log1ru / rho);
        double xi_w = (1.0 + rho * delta_i) / (1.0 + rho * sum_l0u_exp);  // E(xi | zeta, O)
        double denom_increment = ww * firstTerm * secondTerm;
        denom[i] += denom_increment;
        nom[i] += ww * xi_w * exp_mu_new * firstTerm * secondTerm;
      }
    }
  }

  arma::vec result(n);
  for (int i = 0; i < n; i++) {
    result[i] = (denom[i] != 0) ? (nom[i] / denom[i]) : datum::nan;
  }
  return result;
}

//////////////////////////////////////////////////////////////
// Internal precomputation: pieces that DO NOT depend on etaphi.
// Shared across all (q+1)^2 inner Eetaphi evaluations in Hetaphi.
//////////////////////////////////////////////////////////////

struct EetaphiPrecomp {
  std::vector<arma::vec> mu_b;     // per subject, length q
  std::vector<arma::mat> Sigma_b;  // per subject, q x q
  std::vector<arma::vec> Ki_row;   // per subject, length p_len (empty if p_len == 0)
  int n;
  int q;
  int p_len;
};

static EetaphiPrecomp build_precomp(const List& Y, const List& Z,
                                    const List& inv_omega, const List& K,
                                    const arma::mat& D, const List& c_list,
                                    int p_len) {
  int n = Y.size();
  int q = D.n_rows;
  if ((int)D.n_cols != q) stop("D must be square.");

  EetaphiPrecomp pc;
  pc.n = n;
  pc.q = q;
  pc.p_len = p_len;
  pc.mu_b.reserve(n);
  pc.Sigma_b.reserve(n);
  pc.Ki_row.reserve(n);

  for (int i = 0; i < n; i++) {
    arma::vec Yi         = as<arma::vec>(Y[i]);
    arma::mat Zi         = as<arma::mat>(Z[i]);
    arma::mat inv_omegai = as<arma::mat>(inv_omega[i]);
    arma::mat Ki         = as<arma::mat>(K[i]);

    arma::vec prec = inv_omegai.diag();
    if (prec.min() <= 0) stop("inv_omega diagonal must be positive (subject %d).", i + 1);
    arma::mat omega = diagmat(1.0 / prec);

    if ((int)Zi.n_cols != q) {
      stop("Z[%d] has %d cols but q=dim(D)=%d. Z and D dimension mismatch.",
           i + 1, (int)Zi.n_cols, q);
    }

    // V = omega + Z D Z^T   (etaphi-independent)
    arma::mat ZD = Zi * D;
    arma::mat V  = omega + (ZD * Zi.t());

    arma::vec ci       = as<arma::vec>(c_list[i]);
    arma::vec Yi_tilde = Yi - Zi * ci;

    // mu_b = (Z D)^T V^{-1} (Y - Z c) + c
    arma::vec x   = solve(V, Yi_tilde, solve_opts::fast);
    arma::vec mbi = ZD.t() * x + ci;

    // Sigma_b = D - (Z D)^T V^{-1} (Z D)
    arma::mat A   = solve(V, ZD, solve_opts::fast);
    arma::mat Sbi = D - (ZD.t() * A);

    pc.mu_b.push_back(mbi);
    pc.Sigma_b.push_back(Sbi);

    // K[i] cached as length-p_len column vector for dot(pc.Ki_row[i], p)
    arma::vec Ki_rowi;
    if (p_len > 0) {
      if ((int)Ki.n_cols == p_len && (int)Ki.n_rows == 1) {
        Ki_rowi = Ki.row(0).t();        // 1 x p_len  -> p_len
      } else if ((int)Ki.n_rows == p_len && (int)Ki.n_cols == 1) {
        Ki_rowi = Ki.col(0);            // p_len x 1 -> p_len
      } else {
        stop("K[%d] has dim (%d x %d) but phi length is %d; not conformable.",
             i + 1, (int)Ki.n_rows, (int)Ki.n_cols, p_len);
      }
    }
    pc.Ki_row.push_back(Ki_rowi);
  }

  return pc;
}

//////////////////////////////////////////////////////////////
// Internal fast E-step objective that consumes a prebuilt precomp.
// Only the 3 scalars per subject (Kphi, mu_scalar, Sigma_scalar) are
// recomputed per etaphi.
//////////////////////////////////////////////////////////////

static double Eetaphi_fast(const arma::vec& etaphi,
                           const EetaphiPrecomp& pc,
                           const arma::vec& l0i, const List& l0u,
                           const arma::vec& Di, int nK,
                           const arma::vec& w, const arma::vec& v,
                           const List& mu_surv, const List& Sigma2_surv,
                           double rho) {
  int n_etaphi = etaphi.n_elem;
  int q        = pc.q;
  int p_len    = pc.p_len;
  int n        = pc.n;

  if (n_etaphi < q) {
    stop("Length(etaphi) (%d) is smaller than q=dim(D) (%d).", n_etaphi, q);
  }

  arma::vec e = etaphi.subvec(0, q - 1);
  arma::vec p;
  if (p_len > 0) {
    p = etaphi.subvec(q, n_etaphi - 1);
  } else {
    p.set_size(0);
  }

  List mu_surv_new(n);
  List Sigma2_surv_new(n);

  for (int i = 0; i < n; i++) {
    const arma::vec& mbi = pc.mu_b[i];
    const arma::mat& Sbi = pc.Sigma_b[i];

    double Kphi_scalar = 0.0;
    if (p_len > 0) {
      Kphi_scalar = dot(pc.Ki_row[i], p);
    }

    double mu_scalar = Kphi_scalar + dot(e, mbi);
    arma::vec mu1(1); mu1[0] = mu_scalar;
    mu_surv_new[i] = mu1;

    double Sigma_scalar = as_scalar(e.t() * Sbi * e);
    arma::vec s1(1); s1[0] = Sigma_scalar;
    Sigma2_surv_new[i] = s1;
  }

  arma::vec Es_exp = Esurv_exp(w, v, mu_surv, Sigma2_surv,
                               mu_surv_new, Sigma2_surv_new, l0i, l0u, rho);

  arma::vec Es     = Esurv(w, v, mu_surv, Sigma2_surv,
                           mu_surv_new, Sigma2_surv_new, l0i, l0u, rho);

  double result = 0.0;
  for (uword i = 0; i < Di.n_elem; i++) {
    arma::vec l0u_i = l0u[i];
    result += Di[i] * Es[i] - arma::sum(l0u_i) * Es_exp[i];
  }

  return result;
}

static arma::vec Setaphi_core(const arma::vec& etaphi,
                              const EetaphiPrecomp& pc,
                              const arma::vec& l0i, const List& l0u,
                              const arma::vec& Di, int nK,
                              const arma::vec& w, const arma::vec& v,
                              const List& mu_surv, const List& Sigma2_surv,
                              double eps, double rho) {
  int ep_size = etaphi.n_elem;
  arma::vec gradient(ep_size);

  double f0 = Eetaphi_fast(etaphi, pc, l0i, l0u, Di, nK,
                           w, v, mu_surv, Sigma2_surv, rho);

  for (int i = 0; i < ep_size; ++i) {
    arma::vec ep = etaphi;
    double xi = std::max(std::abs(ep[i]), 1.0);
    double h  = xi * eps;

    ep[i] += h;

    double f1 = Eetaphi_fast(ep, pc, l0i, l0u, Di, nK,
                             w, v, mu_surv, Sigma2_surv, rho);

    gradient[i] = (f1 - f0) / h;
  }

  return gradient;
}

//////////////////////////////////////////////////////////////
// Exported wrappers (signatures unchanged; build precomp once then delegate).
//////////////////////////////////////////////////////////////

// [[Rcpp::export]]
double Eetaphi(arma::vec etaphi,
               List Y, List Z, List inv_omega, List K,
               arma::mat D,
               arma::vec l0i, List l0u,
               arma::vec Di, int nK,
               arma::vec w, arma::vec v,
               List mu_surv, List Sigma2_surv,
               List c_list, double rho = 0.0) {

  int q     = D.n_rows;
  int p_len = etaphi.n_elem - q;
  if (p_len < 0) {
    stop("Length(etaphi) (%d) is smaller than q=dim(D) (%d).",
         (int)etaphi.n_elem, q);
  }

  EetaphiPrecomp pc = build_precomp(Y, Z, inv_omega, K, D, c_list, p_len);
  return Eetaphi_fast(etaphi, pc, l0i, l0u, Di, nK, w, v, mu_surv, Sigma2_surv, rho);
}

// [[Rcpp::export]]
arma::vec Setaphi(arma::vec etaphi,
                  List Y, List Z, List inv_omega,
                  List K, arma::mat D,
                  arma::vec l0i, List l0u, arma::vec Di,
                  int nK, arma::vec w, arma::vec v,
                  List mu_surv, List Sigma2_surv, double eps,
                  List c_list, double rho = 0.0) {

  int q     = D.n_rows;
  int p_len = etaphi.n_elem - q;
  if (p_len < 0) {
    stop("Length(etaphi) (%d) is smaller than q=dim(D) (%d).",
         (int)etaphi.n_elem, q);
  }

  EetaphiPrecomp pc = build_precomp(Y, Z, inv_omega, K, D, c_list, p_len);
  return Setaphi_core(etaphi, pc, l0i, l0u, Di, nK, w, v,
                      mu_surv, Sigma2_surv, eps, rho);
}

//////////////////////////////////////////////////////////////
// Per-subject E-step objective: length-n vector of the subject
// contributions  Di[i]*Es[i] - sum(l0u_i)*Es_exp[i]  (NOT summed).
//////////////////////////////////////////////////////////////
static arma::vec Eetaphi_fast_perobs(const arma::vec& etaphi,
                                     const EetaphiPrecomp& pc,
                                     const arma::vec& l0i, const List& l0u,
                                     const arma::vec& Di, int nK,
                                     const arma::vec& w, const arma::vec& v,
                                     const List& mu_surv, const List& Sigma2_surv,
                                     double rho) {
  int n_etaphi = etaphi.n_elem;
  int q        = pc.q;
  int p_len    = pc.p_len;
  int n        = pc.n;
  if (n_etaphi < q) stop("Length(etaphi) (%d) < q (%d).", n_etaphi, q);

  arma::vec e = etaphi.subvec(0, q - 1);
  arma::vec p;
  if (p_len > 0) p = etaphi.subvec(q, n_etaphi - 1); else p.set_size(0);

  List mu_surv_new(n);
  List Sigma2_surv_new(n);
  for (int i = 0; i < n; i++) {
    const arma::vec& mbi = pc.mu_b[i];
    const arma::mat& Sbi = pc.Sigma_b[i];
    double Kphi_scalar = 0.0;
    if (p_len > 0) Kphi_scalar = dot(pc.Ki_row[i], p);
    double mu_scalar = Kphi_scalar + dot(e, mbi);
    arma::vec mu1(1); mu1[0] = mu_scalar; mu_surv_new[i] = mu1;
    double Sigma_scalar = as_scalar(e.t() * Sbi * e);
    arma::vec s1(1); s1[0] = Sigma_scalar; Sigma2_surv_new[i] = s1;
  }

  arma::vec Es_exp = Esurv_exp(w, v, mu_surv, Sigma2_surv,
                               mu_surv_new, Sigma2_surv_new, l0i, l0u, rho);
  arma::vec Es     = Esurv(w, v, mu_surv, Sigma2_surv,
                           mu_surv_new, Sigma2_surv_new, l0i, l0u, rho);

  arma::vec out(n);
  for (int i = 0; i < n; i++) {
    arma::vec l0u_i = l0u[i];
    out[i] = Di[i] * Es[i] - arma::sum(l0u_i) * Es_exp[i];
  }
  return out;
}

//////////////////////////////////////////////////////////////
// Per-subject survival score matrix (n x ep_size): forward-difference
// of the per-subject E-step objective.  colSums() equals Setaphi().
// Used to build the sandwich "meat" sum_i S_i S_i^T = crossprod(U).
//////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat Setaphi_perobs(arma::vec etaphi,
                         List Y, List Z, List inv_omega,
                         List K, arma::mat D,
                         arma::vec l0i, List l0u, arma::vec Di,
                         int nK, arma::vec w, arma::vec v,
                         List mu_surv, List Sigma2_surv, double eps,
                         List c_list, double rho = 0.0) {
  int q     = D.n_rows;
  int p_len = etaphi.n_elem - q;
  if (p_len < 0) {
    stop("Length(etaphi) (%d) is smaller than q=dim(D) (%d).",
         (int)etaphi.n_elem, q);
  }
  EetaphiPrecomp pc = build_precomp(Y, Z, inv_omega, K, D, c_list, p_len);

  int ep_size = etaphi.n_elem;
  int n       = pc.n;
  arma::vec f0 = Eetaphi_fast_perobs(etaphi, pc, l0i, l0u, Di, nK,
                                     w, v, mu_surv, Sigma2_surv, rho);
  arma::mat U(n, ep_size);
  for (int j = 0; j < ep_size; ++j) {
    arma::vec ep = etaphi;
    double xi = std::max(std::abs(ep[j]), 1.0);
    double h  = xi * eps;
    ep[j] += h;
    arma::vec f1 = Eetaphi_fast_perobs(ep, pc, l0i, l0u, Di, nK,
                                       w, v, mu_surv, Sigma2_surv, rho);
    U.col(j) = (f1 - f0) / h;
  }
  return U;
}

// [[Rcpp::export]]
arma::mat Hetaphi(arma::vec etaphi,
                  List Y, List Z, List inv_omega,
                  List K, arma::mat D,
                  arma::vec l0i, List l0u, arma::vec Di,
                  int nK, arma::vec w, arma::vec v,
                  List mu_surv, List Sigma2_surv, double eps,
                  List c_list, double rho = 0.0) {

  int q       = D.n_rows;
  int ep_size = etaphi.n_elem;
  int p_len   = ep_size - q;
  if (p_len < 0) {
    stop("Length(etaphi) (%d) is smaller than q=dim(D) (%d).",
         ep_size, q);
  }

  // Build per-subject precomputation ONCE; reused for (ep_size+1)^2 calls.
  EetaphiPrecomp pc = build_precomp(Y, Z, inv_omega, K, D, c_list, p_len);

  arma::mat hessian(ep_size, ep_size);

  arma::vec f0 = Setaphi_core(etaphi, pc, l0i, l0u, Di, nK, w, v,
                              mu_surv, Sigma2_surv, eps, rho);

  for (int i = 0; i < ep_size; ++i) {
    arma::vec ep = etaphi;
    double xi = std::max(std::abs(ep[i]), 1.0);
    double h  = xi * eps;

    ep[i] += h;

    arma::vec f1 = Setaphi_core(ep, pc, l0i, l0u, Di, nK, w, v,
                                mu_surv, Sigma2_surv, eps, rho);

    hessian.col(i) = (f1 - f0) / h;
  }

  hessian = 0.5 * (hessian + hessian.t());
  return hessian;
}
