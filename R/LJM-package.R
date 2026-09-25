#' @keywords internal
#'
#' @description
#' The local joint model (LJM) for dynamic prediction of survival
#' outcomes from noisy longitudinal biomarkers. LJM uses a local-linear
#' approximation of each biomarker trajectory around a landmark time \eqn{s} and
#' jointly estimates the effect of the recent \emph{level} (intercept) and
#' recent \emph{change} (slope) of the biomarkers on the conditional survival
#' risk, with optional time-varying association functions.
#'
#' @details
#' Typical workflow:
#' \enumerate{
#'   \item \code{\link{LJM_dat}} — build the landmarked longitudinal + survival
#'     dataset at a landmark time \code{s} (with the default complete-case filter).
#'   \item \code{\link{LJM}} — fit the model; supply \code{Bs} for a time-varying
#'     (B-spline) association model.
#'   \item \code{\link[=predict.LJM]{predict}} — conditional survival / risk at a
#'     horizon \code{tau}.
#'   \item \code{\link{AUCdyn}}, \code{\link{PEdyn}} — evaluate time-dependent AUC
#'     and prediction error.
#'   \item \code{\link{confBands.LJM}} — pointwise confidence bands for the
#'     time-varying coefficient functions.
#'   \item \code{\link{select_h_longitudinal}} — cross-validated bandwidth
#'     selection for the local-linear step.
#' }
#' See \code{\link{LJM}} for the model specification and the example dataset
#' \code{\link{pbc2}} for a worked illustration.
#'
#' @useDynLib LJM, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import dplyr
#' @import survival
#' @importFrom splines bs
#' @importFrom nlme lme fixef ranef
#' @importFrom statmod gauss.quad
#' @importFrom rlang :=
#' @importFrom stats as.formula ave complete.cases dnorm model.response predict
#'   qchisq quantile sd var optim optimize pnorm rnorm runif uniroot integrate
#'   cov cor coef vcov na.omit setNames
#' @importFrom utils combn head tail
"_PACKAGE"

# NSE column names used inside dplyr pipelines (silence R CMD check notes)
utils::globalVariables(c(".", "id", "time", "start", "n_eta"))
