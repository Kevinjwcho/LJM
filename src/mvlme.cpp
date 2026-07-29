// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// The first to be used inside a loop (i.e. list indexing done in R)...
// [[Rcpp::export]]
double rcpp_e(const arma::colvec Yik, const arma::mat Xik, const arma::mat Zik, const arma::vec beta, const arma::vec b){
  arma::colvec out = Yik - (Xik * beta + Zik * b);
	return as_scalar(out.t() * out);
}

// The next does the looping internally...
// OLD VERSION - NOT MV EXTENSION!!
/* NumericVector Ee(const Rcpp::List Y, const Rcpp::List X, const Rcpp::List Z,
               const mat& beta, const Rcpp::List b, const int ids, const int K){
	mat M = zeros<mat>(ids, K);
	Rcpp::NumericVector e(K);
	for(int i = 0; i < ids; i++){
	  mat Yi = as<mat>(Y[i]);
	  mat bi = as<mat>(b[i]);
	  Rcpp::List Xi = X[i];
	  Rcpp::List Zi = Z[i];

		for(int k = 0; k < K; k++){
		  colvec Yik = Yi.col(k);
		  mat Xik = as<mat>(Xi[k]);
		  mat Zik = as<mat>(Zi[k]);
		  colvec temp = Yik - (Xik * beta.row(k).t() + Zik * bi.row(k).t());
     	  M(i,k) = as_scalar(temp.t() * temp);
		}
	}
	// Work out column sums
	for(int j = 0; j < K; j++){
		e[j] = sum(M.col(j));
	}

	return e;
} */

// The next does the looping internally...
// [[Rcpp::export]]
NumericVector Ee(const Rcpp::List Y, const Rcpp::List X, const Rcpp::List Z,
               const arma::mat beta, const Rcpp::List b, const Rcpp::List bbT, const int ids, const int K){
  arma::mat M = zeros<mat>(ids, K);
	Rcpp::NumericVector e(K);
	for(int i = 0; i < ids; i++){
	  arma::mat Yi = as<mat>(Y[i]);
	  // Rcpp::Rcout << "Yi size " << Yi.n_rows << ", " << Yi.n_cols << std::endl;
	  arma::mat bi = as<mat>(b[i]);
	  Rcpp::List Xi = X[i];
	  Rcpp::List Zi = Z[i];
	  Rcpp::List bbTi = bbT[i];
	  // Rcpp::Rcout << "Processing ID " << i << std::endl;
		for(int k = 0; k < K; k++){
		  arma::colvec Yik = Yi.col(k);
		  arma::mat Xik = as<mat>(Xi[k]);
		  arma::mat Zik = as<mat>(Zi[k]);
		  arma::mat bbTik = as<mat>(bbTi[k]);

		  // Add debugging output for intermediate variables
		  // Rcpp::Rcout << "Processing K " << k << ": Yik size " << Yik.n_rows << ", " << Yik.n_cols << std::endl;
		  // Rcpp::Rcout << "Xik size " << Xik.n_rows << ", " << Xik.n_cols << std::endl;
		  // Rcpp::Rcout << "Zik size " << Zik.n_rows << ", " << Zik.n_cols << std::endl;
		  // Rcpp::Rcout << "bbTik size " << bbTik.n_rows << ", " << bbTik.n_cols << std::endl;

		  arma::colvec Resid = Yik - Xik * beta.row(k).t();
		  arma::colvec temp = Yik - (Xik * beta.row(k).t() + Zik * bi.row(k).t());
     	  M(i,k) = as_scalar(
			Resid.t() * (Resid - 2.0 * (Zik * bi.row(k).t())) + trace(Zik.t() * Zik * bbTik)
		  );
		}
	}
	// Work out column sums
	for(int j = 0; j < K; j++){
		e[j] = sum(M.col(j));
	}

	return e;
}

// For E[b]
// [[Rcpp::export]]
List Eb(const List Y, const List X, const List Z, const List V, const arma::mat D, const arma::vec beta, const int ids){
	List b(ids);
	for(int i = 0; i < ids; i++){
	  arma::colvec Yi = Y[i];
	  arma::mat Xi = as<mat>(X[i]);
	  arma::mat Zi = as<mat>(Z[i]);
	  arma::mat Vi = as<mat>(V[i]);
	  arma::mat ZD = Zi * D;
	  arma::mat ZDZtV = Zi * D * Zi.t() + Vi;
		b[i] = ZD.t() * inv(ZDZtV) * (Yi - Xi * beta);
	}
	return b;
}

// VarCorr for b
// [[Rcpp::export]]
List covb(const List Z, const List V, const arma::mat Dinv, const int ids){
	List S(ids);
	for(int i = 0; i < ids; i++){
	  arma::mat Zi = as<mat>(Z[i]);
	  arma::mat Vi = as<mat>(V[i]);
		S[i] = inv(Zi.t() * inv(Vi) * Zi + Dinv);
	}
	return S;
}

// E[bbT]
// [[Rcpp::export]]
List EbbT(const List b, const List Sigma, const int ids){
	List bbT(ids);
	for(int i = 0; i < ids;  i++){
	  arma::colvec bi = b[i];
	  arma::mat Sigmai = as<mat>(Sigma[i]);
		bbT[i] = Sigmai + bi * bi.t();
	}
	return bbT;
}

// RHS for beta update
// [[Rcpp::export]]
List betaRHS(const List X, const List Y, const List Z, const List b, const int ids){
	List beta(ids);
	for(int i = 0; i < ids; i++){
	  arma::mat Xi = as<mat>(X[i]);
	  arma::mat Zi = as<mat>(Z[i]);
	  arma::colvec bi = b[i];
	  arma::colvec Yi = Y[i];
		beta[i] = Xi.t() * (Yi - Zi * bi);
	}
	return beta;
}
