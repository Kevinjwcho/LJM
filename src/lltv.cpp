#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// Joint log-likelihood for minimization (-likelihood)
// [[Rcpp::export]]
double ll_tv(arma::vec b, const arma::colvec Y, const arma::mat X, const arma::mat Z, const arma::mat V, const arma::mat D,
          int mi, const arma::rowvec K, const int Delta, const double l0i,
          const arma::colvec l0u,
          const arma::vec beta, const arma::vec phi,
          const arma::mat etaBB,
          const int nK, const int q) {
  try {
    // Debugging dimensions
    // Rcpp::Rcout << "Dimensions of input variables:\n";
    // Rcpp::Rcout << "b: " << b.n_elem << "\nY: " << Y.n_elem << "\nX: " << X.n_rows << "x" << X.n_cols
    //             << "\nZ: " << Z.n_rows << "x" << Z.n_cols << "\nV: " << V.n_rows << "x" << V.n_cols
    //             << "\nD: " << D.n_rows << "x" << D.n_cols << "\n";
    // Rcpp::Rcout << "K: " << K.n_cols << "\nphi: " << phi.n_elem
    //             << "\netaBB: " << etaBB.n_rows << "x" << etaBB.n_cols << "\n";

    double temp = 0.0;
    if (Delta == 1) temp = log(l0i);

    arma::colvec resid = Y - X * beta - Z * b;
    // Rcpp::Rcout << "Residuals calculated. resid: " << resid.n_elem << "\n";

    // Handle edge case for etaBB
    arma::rowvec etaBBVi;  // Declare variable for etaBBVi
    if (etaBB.n_rows == 1) {
      etaBBVi = etaBB.row(0);  // Assign the single row directly
      // Rcpp::Rcout << "etaBB has one row. Using the single row.\n";
    } else if (etaBB.n_rows > 1) {
      etaBBVi = etaBB.row(etaBB.n_rows - 1);  // Use the last row
      // Rcpp::Rcout << "etaBB has multiple rows. Using the last row.\n";
    } else {
      Rcpp::stop("etaBB is empty; cannot proceed.");
    }
    // Rcpp::Rcout << "Last row of etaBB extracted. etaBBVi: " << etaBBVi.n_elem << "\n";

    // Compute Kronecker-like term manually for debugging
    arma::colvec exp_term = kron(exp(K * phi),  exp(etaBB * b));
    // Rcpp::Rcout << "exp(K * phi): " << exp(K * phi) << "\n";
    // Rcpp::Rcout << "exp(etaBB * b): " << exp(etaBB * b).t() << "\n";
    // Rcpp::Rcout << "Element-wise multiplication (exp_term): " << exp_term.n_elem << "\n";

    double result = -1.0*as_scalar(-mi / 2.0 * log(2.0 * M_PI) - 0.5 * log(det(V))
                                - 0.5 * resid.t() * V.i() * resid
                                - nK * log(2.0 * M_PI) - 0.5 * log(det(D)) - 0.5 * b.t() * D.i() * b
                                + temp + Delta * (K * phi + as_scalar(etaBBVi * b))
                                - l0u.t() * exp_term);

                                // Rcpp::Rcout << "Final log-likelihood calculated: " << result << "\n";
                                return result;

  } catch (std::exception& e) {
    Rcpp::Rcerr << "Error in ll function: " << e.what() << "\n";
    return NAN;
  }
}


// [[Rcpp::export]]
arma::mat elementwise_multiply(const arma::mat mat, const arma::vec vec) {
  // Check if sizes are compatible
  if (mat.n_rows != vec.n_elem) {
    Rcpp::stop("The number of rows in the matrix must match the length of the vector.");
  }

  // Perform element-wise multiplication
  return mat.each_col() % vec;
}

// Gradient function of -likelihood w.r.t random effects b_i
// [[Rcpp::export]]
arma::colvec gradll_tv(arma::vec b, const arma::colvec Y, const arma::mat X, const arma::mat Z, const arma::mat V, const arma::mat D,
              int mi, const arma::rowvec K, const int Delta, const double l0i,
              const arma::colvec l0u,
              const arma::vec beta, const arma::vec phi,
              const arma::mat etaBB,
              const int nK, const int q) {
  try {
    // Debugging dimensions
    // Rcpp::Rcout << "Dimensions of input variables for gradient:\n";
    // Rcpp::Rcout << "b: " << b.n_elem << "\nY: " << Y.n_elem << "\nX: " << X.n_rows << "x" << X.n_cols
    //             << "\nZ: " << Z.n_rows << "x" << Z.n_cols << "\nV: " << V.n_rows << "x" << V.n_cols
    //             << "\nD: " << D.n_rows << "x" << D.n_cols << "\n";
    // Rcpp::Rcout << "K: " << K.n_cols << "\nphi: " << phi.n_elem
    //             << "\netaBB: " << etaBB.n_rows << "x" << etaBB.n_cols << "\n";

    arma::colvec resid = Y - X * beta - Z * b;
    // Rcpp::Rcout << "Residuals calculated for gradient. resid: " << resid.n_elem << "\n";

    // Handle edge case for etaBB
    arma::rowvec etaBBVi;  // Declare variable for etaBBVi
    if (etaBB.n_rows == 1) {
      etaBBVi = etaBB.row(0);  // Assign the single row directly
      // Rcpp::Rcout << "etaBB has one row. Using the single row.\n";
    } else if (etaBB.n_rows > 1) {
      etaBBVi = etaBB.row(etaBB.n_rows - 1);  // Use the last row
      // Rcpp::Rcout << "etaBB has multiple rows. Using the last row.\n";
    } else {
      Rcpp::stop("etaBB is empty; cannot proceed.");
    }
    // Rcpp::Rcout << "Last row of etaBB extracted and transposed. etaBBVi: " << etaBBVi.n_elem << "\n";

    arma::vec Delta_vec = Delta * ones<vec>(etaBBVi.n_elem);
    // Rcpp::Rcout << "Delta vector created. Delta_vec: " << Delta_vec.n_elem << "\n";

    // Compute exp_term manually for debugging
    arma::colvec exp_term = kron(exp(K * phi), exp(etaBB * b));
    arma::mat grad_term = elementwise_multiply(etaBB, exp_term);
    arma::rowvec whole_term = l0u.t()* grad_term;
    // Rcpp::Rcout << "exp(K * phi): " << exp(K * phi) << "\n";
    // Rcpp::Rcout << "exp(etaBB * b): " << exp(etaBB * b).t() << "\n";
    // Rcpp::Rcout << "Element-wise multiplication (exp_term): " << exp_term.n_elem << "\n";

    arma::colvec grad = -1.0*(Z.t() * V.i() * resid - D.i() * b
                            + Delta_vec % etaBBVi.t() - whole_term.t());

    // Rcpp::Rcout << "Gradient calculated. grad: " << grad.n_elem << "\n";
    return grad;

  } catch (std::exception& e) {
    Rcpp::Rcerr << "Error in gradll function: " << e.what() << "\n";
    return vec();
  }
}

// Joint likelihood: second derivative w.r.t random effects b_i
// [[Rcpp::export]]
arma::mat sdll_tv(arma::vec b, const arma::mat Z, const arma::mat D, const arma::mat V,
         const arma::rowvec K, const arma::vec l0u, const arma::vec phi,
         const arma::mat etaBB, const int nK) {
  try {
    // Debugging input dimensions
    // Rcpp::Rcout << "Debugging input dimensions:\n";
    // Rcpp::Rcout << "b: " << b.n_elem << "\n";
    // Rcpp::Rcout << "Z: " << Z.n_rows << "x" << Z.n_cols << "\n";
    // Rcpp::Rcout << "D: " << D.n_rows << "x" << D.n_cols << "\n";
    // Rcpp::Rcout << "V: " << V.n_rows << "x" << V.n_cols << "\n";
    // Rcpp::Rcout << "K: " << K.n_elem << "\n";
    // Rcpp::Rcout << "l0u: " << l0u.n_elem << "\n";
    // Rcpp::Rcout << "phi: " << phi.n_elem << "\n";
    // Rcpp::Rcout << "etaBB: " << etaBB.n_rows << "x" << etaBB.n_cols << "\n";

    // Handle edge case for etaBB
    arma::rowvec etaBBVi;  // Declare variable for etaBBVi
    if (etaBB.n_rows == 1) {
      etaBBVi = etaBB.row(0);  // Assign the single row directly
      // Rcpp::Rcout << "etaBB has one row. Using the single row.\n";
    } else if (etaBB.n_rows > 1) {
      etaBBVi = etaBB.row(etaBB.n_rows - 1);  // Use the last row
      // Rcpp::Rcout << "etaBB has multiple rows. Using the last row.\n";
    } else {
      Rcpp::stop("etaBB is empty; cannot proceed.");
    }

    // Rcpp::Rcout << "Last row of etaBB extracted. etaBBVi: " << etaBBVi.n_elem << "\n";

    // Debugging intermediate results
    arma::colvec exp_K_phi = exp(K * phi); // Check this is a scalar
    arma::colvec exp_etaBB_b = exp(etaBB * b); // This should be a vector
    // Rcpp::Rcout << "exp(K * phi): " << exp_K_phi << "\n";
    // Rcpp::Rcout << "exp(etaBB * b): " << exp_etaBB_b.t() << "\n";

    // Compute kernel
    int t = etaBB.n_rows;
    int p = etaBB.n_cols;
    arma::cube outer_products(p, p, t);  // Store outer products for each row
    arma::cube kernel(p, p, t);  // Store outer products for each row
    for (int i = 0; i < t; i++) {
      // Rcpp::Rcout << "Row " << i << " of etaBB: " << etaBB.row(i) << "\n";
      arma::rowvec etaBB_t = etaBB.row(i);         // Extract the i-th row
      outer_products.slice(i) = etaBB_t.t() * etaBB_t;  // Compute the outer product
      kernel.slice(i) = l0u(i) * exp_K_phi(0) * exp_etaBB_b(i) * outer_products.slice(i);
    }
    arma::mat kernel_sum = arma::sum(kernel, 2);  // Sum the kernel matrices
    // Rcpp::Rcout << "Kernel computed: " << kernel << "\n";

    // Compute second derivative matrix
    arma::mat result = -1.0 * (Z.t() * V.i() * Z) - D.i() - kernel_sum;

    // Rcpp::Rcout << "Second derivative matrix computed successfully.\n";
    return result;

  } catch (std::exception& e) {
    Rcpp::Rcerr << "Error in sdll function: " << e.what() << "\n";
    return mat();  // Return an empty matrix in case of failure
  }
}

