// Filename: etaphiRcpp.cpp

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Function to replicate eta row-wise and scale matrices
// [[Rcpp::export]]
List scaleMatricesByEta(const List& BBi, const arma::vec& eta) {
  int n = BBi.size();
  List result(n);
  
  for (int i = 0; i < n; ++i) {
    // Convert each list element to Armadillo matrix
    mat B = as<mat>(BBi[i]);
    
    // Check if the matrix is valid
    if (B.n_rows == 0 || B.n_cols == 0) {
      stop("Element of BBi is not a valid matrix or has invalid dimensions.");
    }
    
    // Create a replicated eta matrix
    mat replicatedEta = repmat(eta.t(), B.n_rows, 1); // Transpose eta for row replication
    
    // Element-wise multiplication
    mat scaledMatrix = replicatedEta % B; // Element-wise multiplication in Armadillo
    result[i] = scaledMatrix;
  }
  
  return result;
}

// Function to create the result matrix using column grouping
// [[Rcpp::export]]
arma::mat createResultMatrix(const arma::mat& A, const arma::vec& etaInd) {
  int nRows = A.n_rows;
  int totalCols = A.n_cols;
  
  if (sum(etaInd) != totalCols) {
    stop("The sum of eta.inds must match the number of columns in A.");
  }
  
  int colStart = 0;
  std::vector<mat> results;
  
  for (double eta : etaInd) { // Use double since arma::vec elements are double
    mat subset = A.cols(colStart, colStart + eta - 1); // Extract subset of columns
    colStart += eta;
    
    // Separate odd and even columns
    vec summedOdds = zeros<vec>(nRows);
    vec summedEvens = zeros<vec>(nRows);
    
    for (uword j = 0; j < subset.n_cols; ++j) {
      if (j % 2 == 0) {
        summedOdds += subset.col(j);
      } else {
        summedEvens += subset.col(j);
      }
    }
    
    // Combine into a two-column matrix
    mat combined(nRows, 2);
    combined.col(0) = summedOdds;
    combined.col(1) = summedEvens;
    results.push_back(combined);
  }
  
  // Combine all results into a single matrix
  int resultCols = results.size() * 2;
  mat resultMatrix(nRows, resultCols);
  
  int colIndex = 0;
  for (const mat& res : results) {
    resultMatrix.cols(colIndex, colIndex + 1) = res;
    colIndex += 2;
  }
  
  return resultMatrix;
}

// Main function to scale and process matrices
// [[Rcpp::export]]
List etaBBs_compute(const List& BBi, const arma::vec& eta, const arma::vec& etaInd) {
  // Scale matrices
  List scaledMatrices = scaleMatricesByEta(BBi, eta);
  
  // Process matrices with createResultMatrix
  List processedMatrices(scaledMatrices.size());
  
  for (int i = 0; i < scaledMatrices.size(); ++i) {
    mat scaled = as<mat>(scaledMatrices[i]);
    processedMatrices[i] = createResultMatrix(scaled, etaInd);
  }
  
  return processedMatrices;
}


// Function to compute Esurv (no exponential in muTerm)
// [[Rcpp::export]]
arma::vec Esurv_tv(arma::vec w, arma::vec v, List mu, List variance, 
                   List mu_new, List variance_new, arma::vec l0i, List l0u) {
  
  int n = mu.size();
  int rho = w.n_elem;
  arma::vec nom(n, fill::zeros);
  arma::vec denom(n, fill::zeros);
  
  for (int l = 0; l < rho; l++) {
    double ww = w[l];
    double vv = v[l];
    
    for (int i = 0; i < n; i++) {
      arma::vec mui = mu[i];
      arma::vec variancei = variance[i];
      arma::vec mu_newi = mu_new[i];
      arma::vec variance_newi = variance_new[i];
      
      double mui_val = mui[mui.n_elem - 1];
      double variancei_val = variancei[variancei.n_elem - 1];
      double mu_newi_val = mu_newi[mu_newi.n_elem - 1];
      double variance_newi_val = variance_newi[variance_newi.n_elem - 1];
      
      double exp_mu_Vi = exp(mui_val + sqrt(variancei_val) * vv);
      double mu_new_Vi= mu_newi_val + sqrt(variance_newi_val) * vv;
      
      double firstTerm = (l0i[i] == 0) ? 1.0 : exp_mu_Vi * l0i[i];
      
      arma::vec exp_mu_t = exp(mui + sqrt(variancei) * vv);
      
      arma::vec l0u_i = l0u[i];
      double sum_l0u_exp = as_scalar(l0u_i.t() * exp_mu_t);
      
      double secondTerm = exp(-sum_l0u_exp); // second term needs integral
      
      double denom_increment = ww * firstTerm * secondTerm;
      denom[i] += denom_increment;
      nom[i] += ww * mu_new_Vi * firstTerm * secondTerm; // we do not need to integral for mu_new_Vi
    }
  }
  
  arma::vec result(n);
  for (int i = 0; i < n; i++) {
    if (denom[i] != 0) {
      result[i] = nom[i] / denom[i];
    } else {
      result[i] = datum::nan;  // Handle division by zero appropriately
    }
  }
  
  return result;
}

// Function to compute Esurv_exp (with exponential in muTerm)
// [[Rcpp::export]]
arma::vec Esurv_exp_tv(arma::vec w, arma::vec v, List mu, List variance, 
                       List mu_new, List variance_new, arma::vec l0i, List l0u) {
  
  int n = mu.size();
  int rho = w.n_elem;
  arma::vec nom(n, fill::zeros);
  arma::vec denom(n, fill::zeros);
  
  for (int l = 0; l < rho; l++) {
    double ww = w[l];
    double vv = v[l];
    
    for (int i = 0; i < n; i++) {
      arma::vec mui = mu[i];
      arma::vec variancei = variance[i];
      arma::vec mu_newi = mu_new[i];
      arma::vec variance_newi = variance_new[i];
      
      double mui_val = mui[mui.n_elem - 1];
      double variancei_val = variancei[variancei.n_elem - 1];
      double mu_newi_val = mu_newi[mu_newi.n_elem - 1];
      double variance_newi_val = variance_newi[variance_newi.n_elem - 1];
      
      double exp_mu_Vi = exp(mui_val + sqrt(variancei_val) * vv);
      double mu_new_Vi= mu_newi_val + sqrt(variance_newi_val) * vv;
      
      double firstTerm = (l0i[i] == 0) ? 1.0 : exp_mu_Vi * l0i[i];
      
      arma::vec exp_mu_t = exp(mui + sqrt(variancei) * vv);
      arma::vec exp_mu_new_t = exp(mu_newi + sqrt(variance_newi) * vv);
      
      arma::vec l0u_i = l0u[i];
      double sum_l0u_exp = as_scalar(l0u_i.t() * exp_mu_t);
      double sum_l0u_exp_new = as_scalar(l0u_i.t() * exp_mu_new_t);
      
      double secondTerm = exp(-sum_l0u_exp);
      
      double denom_increment = ww * firstTerm * secondTerm;
      denom[i] += denom_increment;
      nom[i] += ww * sum_l0u_exp_new * firstTerm * secondTerm;
    }
  }
  
  arma::vec result(n);
  for (int i = 0; i < n; i++) {
    if (denom[i] != 0) {
      result[i] = nom[i] / denom[i];
    } else {
      result[i] = datum::nan;  // Handle division by zero appropriately
    }
  }
  
  return result;
}

// [[Rcpp::export]]
List Esurv_exp_t(arma::vec w, arma::vec v, List mu, List variance, 
                 List mu_new, List variance_new, arma::vec l0i, List l0u) {
  
  int n = mu.size();
  int rho = w.n_elem;
  List nom(n);
  arma::vec denom(n, fill::zeros);
  
  for (int l = 0; l < rho; l++) {
    double ww = w[l];
    double vv = v[l];
    
    for (int i = 0; i < n; i++) {
      arma::vec mui = mu[i];
      arma::vec variancei = variance[i];
      arma::vec mu_newi = mu_new[i];
      arma::vec variance_newi = variance_new[i];
      
      double mui_val = mui[mui.n_elem - 1];
      double variancei_val = variancei[variancei.n_elem - 1];
      double mu_newi_val = mu_newi[mu_newi.n_elem - 1];
      double variance_newi_val = variance_newi[variance_newi.n_elem - 1];
      
      double exp_mu_Vi = exp(mui_val + sqrt(variancei_val) * vv);
      double mu_new_Vi= mu_newi_val + sqrt(variance_newi_val) * vv;
      
      double firstTerm = (l0i[i] == 0) ? 1.0 : exp_mu_Vi * l0i[i];
      
      arma::vec exp_mu_t = exp(mui + sqrt(variancei) * vv);
      arma::vec exp_mu_new_t = exp(mu_newi + sqrt(variance_newi) * vv);
      
      arma::vec l0u_i = l0u[i];
      double sum_l0u_exp = as_scalar(l0u_i.t() * exp_mu_t);
      // double sum_l0u_exp_new = as_scalar(l0u_i.t() * exp_mu_new_t);
      
      double secondTerm = exp(-sum_l0u_exp);
      arma::vec nom_sub(exp_mu_new_t.n_elem, fill::zeros);
      double denom_increment = ww * firstTerm * secondTerm;
      denom[i] += denom_increment;
      // Initialize nom[i] if not already set
      if (nom[i] == R_NilValue) {
        nom[i] = exp_mu_new_t * denom_increment;  // Initialize
      } else {
        nom[i] = as<arma::vec>(nom[i]) + exp_mu_new_t * denom_increment;
      }
    }
  }
  
  List result(n);
  for (int i = 0; i < n; i++) {
    if (denom[i] != 0) {
      double denom_i = denom[i];
      arma::vec nom_i = nom[i];
      result[i] = nom_i / denom_i;
    } else {
      result[i] = datum::nan;  // Handle division by zero appropriately
    }
  }
  
  return result;
}

// Function to compute the expected log-likelihood
// [[Rcpp::export]]
double Eetaphi_t(arma::vec etaphi, List X, List Y, List Z, List inv_omega, List K, 
                 arma::mat D, arma::vec beta, arma::vec l0i, List l0u, 
                 arma::vec Di, int nK, arma::vec w, arma::vec v, 
                 List mu_surv, List Sigma2_surv, const List& BBi, const arma::vec& etaInd) {
  
  int n_etaphi = etaphi.n_elem;
  int n_eta = sum(etaInd);
  arma::vec e = etaphi.subvec(0, n_eta - 1);
  arma::vec p = etaphi.subvec(n_eta, n_etaphi - 1);
  
  List eBBs = etaBBs_compute(BBi, e, etaInd);
  
  List mu_surv_new(X.size());
  List Sigma2_surv_new(X.size());
  
  for (int i = 0; i < X.size(); i++) {
    arma::mat Xi = as<arma::mat>(X[i]);
    arma::vec Yi = as<arma::vec>(Y[i]);
    arma::mat Zi = as<arma::mat>(Z[i]);
    arma::mat inv_omegai = as<arma::mat>(inv_omega[i]);
    arma::mat Ki = as<arma::mat>(K[i]);
    arma::mat eBBi = as<arma::mat>(eBBs[i]);
    
    // Compute omega and mu_b
    arma::mat omega = inv(inv_omegai);
    arma::mat V = omega + Zi * D * Zi.t();
    arma::vec mu_b = D * Zi.t() * solve(V, Yi - Xi * beta);
    
    // Compute mu_surv_i
    arma::vec mu_surv_i = ones<vec>(eBBi.n_rows)*(Ki * p) + eBBi * mu_b;
    mu_surv_new[i] = mu_surv_i;
    
    // Compute Sigma_surv_i
    arma::mat Sigma_b = D - D * Zi.t() * solve(V, Zi * D);
    arma::mat Sigma_surv_i = eBBi * Sigma_b * eBBi.t();
    // Extract diagonal of Sigma_surv and convert to arma::vec
    arma::vec diag_sigma_surv = arma::vec(Sigma_surv_i.diag()); // Convert diagview to vec
    Sigma2_surv_new[i] = wrap(diag_sigma_surv);
  }
  
  // Compute expectations using Esurv and Esurv_exp
  arma::vec Es_exp = Esurv_exp_tv(w, v, mu_surv, Sigma2_surv, 
                                  mu_surv_new, Sigma2_surv_new, l0i, l0u);
  
  arma::vec Es = Esurv_tv(w, v, mu_surv, Sigma2_surv, 
                          mu_surv_new, Sigma2_surv_new, l0i, l0u);
  
  // Compute the expected log-likelihood
  double result = 0.0;
  for (unsigned int i = 0; i < Di.n_elem; i++) {
    arma::vec l0u_i = l0u[i];
    result += Di[i] * Es[i] - Es_exp[i];
  }
  
  return result;
}

// Function to compute the gradient of Eetaphi
// [[Rcpp::export]]
arma::vec Setaphi_t(arma::vec etaphi, 
                    List X, List Y, List Z, List inv_omega, 
                    List K, arma::mat D, arma::vec beta, 
                    arma::vec l0i, List l0u, arma::vec Di, 
                    int nK, arma::vec w, arma::vec v, 
                    List mu_surv, List Sigma2_surv,
                    const List& BBi, const arma::vec& etaInd, double eps) {
  
  int ep_size = etaphi.n_elem;
  arma::vec gradient(ep_size);
  double f0 = Eetaphi_t(etaphi, X, Y, Z, inv_omega, K, D, beta, l0i, l0u, Di, nK, w, v, mu_surv, Sigma2_surv, BBi, etaInd);
  
  for (int i = 0; i < ep_size; ++i) {
    arma::vec ep = etaphi;
    double xi = std::max(std::abs(ep[i]), 1.0);
    double h = xi * eps;
    ep[i] += h;
    double f1 = Eetaphi_t(ep, X, Y, Z, inv_omega, K, D, beta, l0i, l0u, Di, nK, w, v, mu_surv, Sigma2_surv, BBi, etaInd);
    gradient[i] = (f1 - f0) / h;
  }
  
  return gradient;
}

// Function to compute the Hessian of Eetaphi
// [[Rcpp::export]]
arma::mat Hetaphi_t(arma::vec etaphi, 
                    List X, List Y, List Z, List inv_omega, 
                    List K, arma::mat D, arma::vec beta, 
                    arma::vec l0i, List l0u, arma::vec Di, 
                    int nK, arma::vec w, arma::vec v, 
                    List mu_surv, List Sigma2_surv,
                    const List& BBi, const arma::vec& etaInd, double eps) {
  
  int ep_size = etaphi.n_elem;
  arma::mat hessian(ep_size, ep_size);
  arma::vec f0 = Setaphi_t(etaphi, X, Y, Z, inv_omega, K, D, beta, 
                           l0i, l0u, Di, nK, w, v, mu_surv, Sigma2_surv,
                           BBi, etaInd, eps);
  
  for (int i = 0; i < ep_size; ++i) {
    arma::vec ep = etaphi;
    double xi = std::max(std::abs(ep[i]), 1.0);
    double h = xi * eps;
    ep[i] += h;
    arma::vec f1 = Setaphi_t(ep, X, Y, Z, inv_omega, K, D, beta, 
                             l0i, l0u, Di, nK, w, v, mu_surv, Sigma2_surv,
                             BBi, etaInd, eps);
    hessian.col(i) = (f1 - f0) / h;
  }
  
  // Ensure the Hessian is symmetric
  hessian = 0.5 * (hessian + hessian.t());
  
  return hessian;
}
