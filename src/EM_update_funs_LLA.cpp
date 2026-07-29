#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//////////////////////////////////////////////////////////////////
// Helpers
//////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
arma::uvec subtract_one_from_uvec(arma::uvec vec) {
  for (arma::uword i = 0; i < vec.size(); i++) {
    if (vec[i] > 0) vec[i] -= 1;  // R (1-based) -> C++ (0-based)
  }
  return vec;
}

//////////////////////////////////////////////////////////////////
// 1) Kernel-weighted inv_omega (precision) for stacked longitudinal vectors
//
// Assumptions:
// - mi_list[[i]]: IntegerVector length nK, counts per marker k
// - var_e: length nK, marker-specific residual variances
// - w_stack_list[[i]]: numeric vector length sum(mi_list[[i]]),
//   aligned with the SAME stacking order used for Y/Z in your code:
//   marker 1 repeated mi[1] times, then marker 2 repeated mi[2] times, ...
//////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List calc_inv_omega_kernel(List mi_list, arma::vec var_e, List w_stack_list) {
  int n = mi_list.size();
  List inv_omega_list(n);
  
  if ((int)var_e.n_elem <= 0) stop("var_e must be non-empty.");
  
  for (int i = 0; i < n; i++) {
    IntegerVector mi = mi_list[i];
    int nK = mi.size();
    
    if ((int)var_e.n_elem != nK) {
      stop("var_e length (%d) must match mi length (%d) for subject %d.",
           (int)var_e.n_elem, nK, i + 1);
    }
    
    int total_mi = sum(mi);
    arma::vec w = as<arma::vec>(w_stack_list[i]);
    
    if ((int)w.n_elem != total_mi) {
      stop("w_stack_list[%d] length (%d) must equal sum(mi_list[%d]) (%d).",
           i + 1, (int)w.n_elem, i + 1, total_mi);
    }
    
    arma::vec diag_prec(total_mi);
    int idx = 0;
    
    // Stacking order matches your original calc_inv_omega:
    // marker k value repeated mi[k] times.
    for (int k = 0; k < nK; k++) {
      if (var_e[k] <= 0) stop("var_e[%d] must be > 0.", k + 1);
      double prec_k = 1.0 / var_e[k];
      
      for (int r = 0; r < mi[k]; r++) {
        diag_prec[idx] = w[idx] * prec_k; // kernel-weighted precision
        idx++;
      }
    }
    
    inv_omega_list[i] = arma::diagmat(diag_prec);
  }
  
  return inv_omega_list;
}

//////////////////////////////////////////////////////////////////
// 2) Sigma.longK = diag(Z * S * Z.t()) for each subject and marker
//
// Inputs:
// - Z_list: list of length n (subjects), each element is a list length nK of matrices Z_{ik}
// - S_list: list of length n (subjects), each element is a list length nK of matrices S_{ik}
// Output:
// - Sigma_longK: list length n, each element is a list length nK of vectors diag(Z S Z^T)
//////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List calc_Sigma_longK(List Z_list, List S_list, int nK) {
  int n = S_list.size();
  List Sigma_longK(n);
  
  for (int j = 0; j < n; j++) {
    List Z_i = Z_list[j];
    List S_i = S_list[j];
    
    if (Z_i.size() != nK || S_i.size() != nK) {
      stop("Subject %d: Z_list and S_list must each have length nK=%d.", j + 1, nK);
    }
    
    List out(nK);
    
    for (int k = 0; k < nK; k++) {
      arma::mat Zk = Z_i[k];
      arma::mat Sk = S_i[k];
      
      out[k] = arma::diagvec(Zk * Sk * Zk.t()); // vector of diag elements
    }
    
    Sigma_longK[j] = out;
  }
  
  return Sigma_longK;
}

//////////////////////////////////////////////////////////////////
// 3) mu.longK without X/beta: mu_{ik} = Y_{ik} - Z_{ik} * b_i(subset)
//
// Inputs:
// - Y_list: list length n, each is matrix with nK columns (marker-wise), column k = Y_{ik}
// - Z_list: list length n, each is list length nK of matrices Z_{ik}
// - b_list: list length n, each is vector b_i
// - b_inds: list length nK, each is uvec of indices (1-based from R) selecting elements of b_i for marker k
//////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List calc_mu_longK_noX(List Y_list, List Z_list, List b_list, List b_inds, int nK) {
  int n = Y_list.size();
  List mu_longK(n);
  
  if (b_inds.size() != nK) stop("b_inds must have length nK.");
  
  for (int i = 0; i < n; i++) {
    arma::mat Y_i = Y_list[i];
    List Z_i = Z_list[i];
    arma::vec b_i = b_list[i];
    
    if ((int)Y_i.n_cols != nK) {
      stop("Subject %d: Y_i must have nK=%d columns.", i + 1, nK);
    }
    if (Z_i.size() != nK) {
      stop("Subject %d: Z_list element must have length nK=%d.", i + 1, nK);
    }
    
    List out(nK);
    
    for (int k = 0; k < nK; k++) {
      arma::vec Yk = Y_i.col(k);
      arma::mat Zk = Z_i[k];
      
      arma::uvec b_ind = b_inds[k];
      arma::uvec b0 = subtract_one_from_uvec(b_ind);
      
      if (b0.max() >= b_i.n_elem) {
        stop("Subject %d marker %d: b_inds out of bounds for b_i.", i + 1, k + 1);
      }
      
      out[k] = Yk - Zk * b_i.elem(b0);
    }
    
    mu_longK[i] = out;
  }
  
  return mu_longK;
}

//////////////////////////////////////////////////////////////////
// 4) Kernel-weighted Ee for each subject and marker
//
// Ee_{ik} = sum_r w_{ikr} * ( mu_{ikr}^2 + SigmaDiag_{ikr} )
//
// Inputs:
// - mu_list: list length n, each is list length nK of vectors mu_{ik}
// - Sigma_longK: list length n, each is list length nK of vectors diag(Z S Z^T)
// - w_longK: list length n, each is list length nK of weight vectors aligned with mu_{ik}
// Output:
// - Ee: list length n, each is rowvec length nK
//////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List calc_Ee_kernel(List mu_list, List Sigma_longK, List w_longK, int nK) {
  int n = mu_list.size();
  List Ee(n);
  
  for (int i = 0; i < n; i++) {
    List mu_i = mu_list[i];
    List Sig_i = Sigma_longK[i];
    List w_i  = w_longK[i];
    
    if (mu_i.size() != nK || Sig_i.size() != nK || w_i.size() != nK) {
      stop("Subject %d: mu_list, Sigma_longK, w_longK must each have length nK=%d.", i + 1, nK);
    }
    
    arma::rowvec out(nK);
    
    for (int k = 0; k < nK; k++) {
      arma::vec muik   = as<arma::vec>(mu_i[k]);
      arma::vec sigd   = as<arma::vec>(Sig_i[k]); // diag vector
      arma::vec wik    = as<arma::vec>(w_i[k]);
      
      if (muik.n_elem != sigd.n_elem || muik.n_elem != wik.n_elem) {
        stop("Subject %d marker %d: length mismatch among mu, SigmaDiag, weights.",
             i + 1, k + 1);
      }
      
      // sum w * (mu^2 + SigmaDiag)
      out[k] = arma::as_scalar( arma::sum( wik % (muik % muik + sigd) ) );
    }
    
    Ee[i] = out;
  }
  
  return Ee;
}