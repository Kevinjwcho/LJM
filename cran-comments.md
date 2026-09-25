## Submission summary

This is the first CRAN submission of `LJM`, an R package for dynamic survival
prediction from longitudinal biomarkers using the local joint model (compiled C++ via Rcpp / RcppArmadillo / RcppEigen).

## Test environments

* local macOS, R 4.5.1
* (please also run: win-builder devel/release, R-hub, GitHub Actions
  ubuntu/windows/macOS before submission)

## R CMD check results

`R CMD check --as-cran` gives 0 ERRORs and 0 WARNINGs.

Remaining NOTEs (with justification):

* "no visible global function definition" for a few internal Rcpp routines
  (e.g. `calc_Ee`, `calc_inv_omega`): these are C++ functions exported through
  `RcppExports` and are defined; the note is the usual static-analysis
  false positive for Rcpp-registered functions called from R.
* "no visible binding for global variable" for column names used inside dplyr
  pipelines (`.`, `id`, `time`, ...); these are registered via
  `utils::globalVariables()`.

If an installation warning of the form `unknown warning group
'-Wfixed-enum-extension'` appears, it originates from R's own
`R_ext/Boolean.h` header under a newer clang than the one R was built with; it
is benign and does not occur on the CRAN build machines.

## Downstream dependencies

There are no reverse dependencies (first submission).
