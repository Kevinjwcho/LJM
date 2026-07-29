#' Fit an LLA-based Jointly Estimated Landmarking (JEL) model
#'
#' Fits a landmarked joint model using Local Linear Approximation (LLA) for the
#' longitudinal component and a Cox model for the survival component, following
#' a `prep -> fit -> return` workflow producing a `"JEL"` object.
#'
#' This function assumes the following helpers exist in the package namespace:
#' \itemize{
#'   \item \code{prep_LLA_landmark()}
#'   \item \code{make_cox_formula_from_LLA()}
#'   \item \code{InitVal_LLAJEL()}
#'   \item \code{RefinedfastEM_LLA()}
#' }
#'
#' @param train_dataset A list containing at least \code{LMM_dat}, \code{Surv_dat},
#'   and \code{var_list}. \code{Surv_dat} must already be in (start, stop, event)
#'   form appropriate for \code{Surv(start, stop, event)}.
#' @param y_vars Character vector of longitudinal outcome names to be modeled by LLA.
#' @param s Landmark time (scalar).
#' @param h Bandwidth for kernel weighting in the LLA landmarking step (scalar).
#' @param base_terms Character vector (or formula components) for baseline Cox terms
#'   to include in the survival submodel alongside LLA-derived terms.
#' @param ker Kernel type for landmark weighting. Default is \code{"gaussian"}.
#' @param model Longitudinal model type: \code{"Fixed"} for constant landmark
#'   associations or \code{"TimeVar"} for time-varying (B-spline) associations.
#'   Default is \code{"Fixed"}.
#' @param Bs Optional list of B-spline specifications (one per longitudinal
#'   outcome) used only when \code{model = "TimeVar"}; each element is a list with
#'   entries \code{df}, \code{degree}, \code{knots}, \code{Bknots} (or \code{NULL}
#'   to leave a marker time-invariant). Default is \code{NULL} (all fixed).
#' @param gh.nodes Number of Gauss-Hermite nodes used by the EM step. Default is 5.
#' @param max.iter Maximum number of EM iterations. Default is 200.
#' @param tol Convergence tolerance. Default is 0.01.
#' @param diff.type Convergence criterion type passed to the EM routine. Default is \code{"abs.rel"}.
#' @param collect.hist Logical; whether to collect iteration history. Default is TRUE.
#' @param Vcov Logical; whether to compute/post-process variance-covariance. Default is TRUE.
#' @param verbose Logical; verbosity for the EM routine. Default is FALSE.
#' @param h_init Optional bandwidth used only to construct the initial local-linear
#'   BLUPs of the random effects, decoupled from the fitting bandwidth \code{h}.
#'   If \code{NULL} (default) the fitting bandwidth \code{h} is reused for
#'   initialization.
#'
#' @details
#' Internally, this function:
#' \enumerate{
#'   \item Builds landmarked LLA objects via \code{prep_LLA_landmark()}.
#'   \item Constructs and fits a Cox model on \code{prep$surv_s2} using
#'     \code{make_cox_formula_from_LLA()} and \code{survival::coxph()}.
#'   \item Builds indexing/design objects for the EM step.
#'   \item Initializes survival parameters via \code{InitVal_LLAJEL()} and runs
#'     the EM algorithm via \code{RefinedfastEM_LLA()}.
#'   \item Returns an object of class \code{"JEL"} with estimates and optional diagnostics.
#' }
#'
#' \strong{Time-varying associations.} When \code{Bs} is supplied (equivalently
#' \code{model = "TimeVar"}), the constant association vector is replaced by a
#' coefficient function \eqn{\eta(u)} expanded in a B-spline basis,
#' \eqn{\eta(u) = \sum_{m=1}^{M} \eta_m B_m(u)}, and the conditional cumulative
#' hazard at landmark \eqn{s} over a horizon \eqn{\tau} is
#' \deqn{\Lambda(\tau \mid s) = G\left[ \int_0^\tau \lambda_0(u)
#'   \exp\left\{ \mathbf{W}_i^\top \phi + \eta(u)^\top \mathbf{b}_i(s) \right\}
#'   du \right],}
#' with transformation \eqn{G(x) = \log(1+\rho x)/\rho} (\eqn{\rho = 0}: Cox
#' proportional hazards; \eqn{\rho = 1}: proportional odds), baseline hazard
#' \eqn{\lambda_0}, and local intercept/slope random effects \eqn{\mathbf{b}_i(s)}
#' at the landmark. Unlike the constant-coefficient case the cumulative hazard
#' does not factor as \eqn{\Lambda_0(\tau)\exp(\cdot)} and is evaluated
#' numerically. Pointwise confidence bands for the intercept function
#' \eqn{\eta_0(u)} and the slope function \eqn{\eta_1(u)} are available via
#' \code{\link{confBands.JEL}}.
#'
#' @return
#' An object of class \code{"JEL"} (a list) containing at least:
#' \itemize{
#'   \item \code{coefficients}: estimated fixed effects / parameters (from EM output).
#'   \item \code{Vcov}: estimated variance-covariance matrix (if computed).
#'   \item \code{est.bi}: estimated random effects (if returned by EM).
#'   \item \code{convergence}: \code{"success"} or \code{"failure"} based on iteration history.
#'   \item \code{history}: EM iteration history (if collected).
#'   \item \code{n}, \code{d}: subject count and event indicators.
#'   \item \code{fitCOX}: fitted Cox model object (for debugging/inspection).
#'   \item \code{prep}: preprocessing output (for debugging/inspection).
#'   \item \code{dataset}: echoes key inputs and settings.
#' }
#'
#' @examples
#' \donttest{
#' data("pbc2", package = "JEL")
#' d <- pbc2; d$id <- as.numeric(d$id)
#' d$Y.1 <- log(d$serBilir); d$Y.2 <- d$albumin
#' vl <- list(id = "id", time = "year", EvTime = "years", event = "status2")
#'
#' ## 1) landmarked data at s = 5 (complete-case filter applied by default)
#' td <- JEL_dat(d, s = 5, var_list = vl, h = 4, y_vars = c("Y.1", "Y.2"))
#'
#' ## 2) fit (Epanechnikov kernel, fixed bandwidth h = 4)
#' fit <- JEL(td, y_vars = c("Y.1", "Y.2"), s = 5, h = 4,
#'            base_terms = "drug", ker = "epanechnikov")
#' fit$coefficients$eta
#'
#' ## 3) conditional risk at horizon tau = 2, then evaluate
#' pr  <- predict(fit, testdat = d, tau = 2)
#' ev  <- list(id = "id", EvTime = "years", event = "status2")
#' AUCdyn(pr, data = d, landmarks = 5, tau = 2, var_list = ev)
#'
#' ## time-varying association (supply a B-spline basis) + confidence bands
#' spl <- list(df = 3, degree = 1, knots = NULL, Bknots = NULL)
#' fit_tv <- JEL(td, y_vars = c("Y.1", "Y.2"), s = 5, h = 4, base_terms = "drug",
#'               ker = "epanechnikov", Bs = list(spl, list(NULL)))
#' cb <- confBands.JEL(fit_tv, K = 1)
#' }
#'
#' @seealso \code{\link{JEL_dat}} for building the input;
#'   \code{\link[=predict.JEL]{predict}} for prediction;
#'   \code{\link{AUCdyn}}, \code{\link{PEdyn}} for evaluation;
#'   \code{\link{confBands.JEL}} for time-varying confidence bands;
#'   \code{\link{select_h_longitudinal}} for bandwidth selection.
#'
#' @export
#'
#' @importFrom survival coxph
#' @importFrom stats model.frame model.matrix formula terms
JEL <- function(train_dataset,
                   y_vars,
                   s, h,
                   base_terms,
                   ker = "gaussian",
                   model = "Fixed",
                   Bs = NULL,
                   gh.nodes = 5,
                   max.iter = 200,
                   tol = 0.01,
                   diff.type = "abs.rel",
                   collect.hist = TRUE,
                   Vcov = TRUE,
                   verbose = FALSE,
                   h_init = NULL) {

  # Determine model type from Bs
  if (!is.null(Bs)) {
    if (length(Bs) != length(y_vars)) {
      stop("Bs must be a list of length equal to the number of y_vars.")
    }
    model <- "TimeVar"
  }
  
  # 1) prep
  prep <- prep_LLA_landmark(
    LMM_dat  = train_dataset$LMM_dat,
    Surv_dat = train_dataset$Surv_dat,
    y_vars   = y_vars,
    s       = s,
    h        = h,
    ker      = ker,
    var_list = train_dataset$var_list,
    h_init   = h_init
  )

  surv_s2 <- prep$surv_s2
  
  # 2) Cox fit on surv_s2 (time-dependent start/stop/event assumed already in Surv_dat)
  cox_form <- make_cox_formula_from_LLA(
    LLA_list     = prep$LLA_list_s,
    base_terms   = base_terms
  )
  
  missing_vars <- setdiff(all.vars(cox_form), names(surv_s2))
  if (length(missing_vars) > 0) {
    stop(paste0("Missing variables in surv_s2: ", paste(missing_vars, collapse = ", ")))
  }
  
  fitCOX <- coxph(
    formula = cox_form,
    data    = surv_s2,
    x = TRUE, model = TRUE,
    id = id, cluster = id
  )
  
  # 3) Build survival indexing objects from fitCOX (like your original code)
  if (is.null(fitCOX$y) || ncol(fitCOX$y) != 3) {
    stop("fitCOX must be specified with Surv(start, stop, event).")
  }
  start <- fitCOX$y[, 1]
  stop  <- fitCOX$y[, 2]
  event <- fitCOX$y[, 3]
  
  ID1_surv <- as.vector(surv_s2$id)
  uniqueID_surv <- !duplicated(ID1_surv)
  tempID_surv <- which(uniqueID_surv)
  tempID_surv <- c(tempID_surv, length(ID1_surv) + 1)
  ni_surv <- diff(tempID_surv)
  ID_surv <- rep(seq_len(sum(uniqueID_surv)), times = ni_surv)
  
  Time <- stop[cumsum(ni_surv)]
  d <- event[cumsum(ni_surv)]
  n_subject <- length(Time)
  
  U <- sort(unique(Time[d == 1]))
  tempU <- lapply(Time, function(t) U[t >= U])
  times <- unlist(tempU)
  nk <- sapply(tempU, length)
  M <- sum(nk)
  
  Indcs <- list(
    Index  = rep(seq_len(n_subject), nk),
    Index0 = match(Time, U),
    Index1 = unlist(lapply(nk[nk != 0], seq, from = 1)),
    Index2 = colSums(d * outer(Time, U, "=="))
  )
  
  # Wtime and Wtime2 from Cox design matrix
  TermsSurv <- fitCOX$terms
  formSurv  <- formula(fitCOX)
  ncw <- ncol(as.matrix(fitCOX$x))
  
  mfSurv <- model.frame(TermsSurv, surv_s2)[cumsum(ni_surv), , drop = FALSE]
  if (ncw > 0) {
    Wtime <- model.matrix(formSurv, mfSurv)
    if (attr(TermsSurv, "intercept")) Wtime <- Wtime[, -1, drop = FALSE]
  } else {
    Wtime <- matrix(nrow = n_subject, ncol = 0)
  }
  
  mfSurv2 <- mfSurv[Indcs$Index, , drop = FALSE]
  if (ncw > 0) {
    Wtime2 <- model.matrix(formSurv, mfSurv2)
    if (attr(TermsSurv, "intercept")) Wtime2 <- Wtime2[, -1, drop = FALSE]
  } else {
    Wtime2 <- matrix(nrow = M, ncol = 0)
  }
  
  # Build subject-split lists (Z.st, Y.st, Wker.st) from prep
  Z.st <- lapply(seq_along(y_vars), function(k) {
    split(prep$Z[[k]], prep$ID[[k]]) |>
      lapply(function(v) matrix(v, ncol = ncol(prep$Z[[k]])))
  })
  Y.st <- lapply(seq_along(y_vars), function(k) split(prep$Y[[k]], prep$ID[[k]]))
  Wker.st <- lapply(seq_along(y_vars), function(k) split(prep$W_ker[[k]], prep$ID[[k]]))

  # 4) Init values + EM fit (branch on model type)
  if (model == "TimeVar") {
    # ---- Time-Varying path ----
    surv.init <- InitVal_LLAJEL_tv(
      bBLUP      = prep$bLLA,
      n_subject  = n_subject,
      n_LME      = length(y_vars),
      Indcs      = Indcs,
      start      = start, stop = stop, event = event,
      W          = Wtime, ncw = ncw,
      Wtime2     = Wtime2,
      ID_surv    = ID_surv,
      Bs         = Bs
    )

    B_list <- surv.init$B_list
    eta_n  <- surv.init$eta_n

    theta <- RefinedfastEM_LLA_tv(
      data       = train_dataset$LMM_dat,
      n_LLA      = length(y_vars),
      ID_surv    = ID_surv,
      ni_surv    = ni_surv,
      Z.st       = Z.st,
      Y.st       = Y.st,
      Wker.st    = Wker.st,
      Wtime      = Wtime,
      Ysigma2    = prep$Ysigma2,
      Bsigma     = prep$Bsigma,
      bLLA       = prep$bLLA,
      cLLA       = prep$cLLA,
      n_subject  = n_subject,
      surv.init  = surv.init,
      ni         = prep$ni,
      ncw        = ncw,
      B_list     = B_list,
      eta_n      = eta_n,
      gh.nodes   = gh.nodes,
      collect.hist = collect.hist,
      max.iter   = max.iter,
      tol        = tol,
      diff.type  = diff.type,
      post.process = Vcov,
      verbose    = verbose
    )
  } else {
    # ---- Fixed (original) path ----
    surv.init <- InitVal_LLAJEL(
      bBLUP      = prep$bLLA,
      n_subject  = n_subject,
      n_LME      = length(y_vars),
      Indcs      = Indcs,
      start      = start, stop = stop, event = event,
      W          = Wtime, ncw = ncw,
      Wtime2     = Wtime2,
      ID_surv    = ID_surv
    )

    theta <- RefinedfastEM_LLA(
      data       = train_dataset$LMM_dat,
      n_LLA      = length(y_vars),
      ID_surv    = ID_surv,
      ni_surv    = ni_surv,
      Z.st       = Z.st,
      Y.st       = Y.st,
      Wker.st    = Wker.st,
      Wtime      = Wtime,
      Ysigma2    = prep$Ysigma2,
      Bsigma     = prep$Bsigma,
      bLLA       = prep$bLLA,
      cLLA       = prep$cLLA,
      n_subject  = n_subject,
      surv.init  = surv.init,
      ni         = prep$ni,
      ncw        = ncw,
      gh.nodes   = gh.nodes,
      collect.hist = collect.hist,
      max.iter   = max.iter,
      tol        = tol,
      diff.type  = diff.type,
      post.process = Vcov,
      verbose    = verbose
    )
  }

  # 5) Return a JEL-like object + attach prep/fitCOX for debugging
  out <- list(
    coefficients = theta$coeffs,
    Vcov         = theta$Vcov,
    Hessian      = theta$Hessian,   # diagnostic passthrough (raw PRES Hessian, tv)
    est.bi       = theta$REs,
    convergence  = if (!is.null(theta$history) && nrow(theta$history) < max.iter) "success" else "failure",
    history      = theta$history,
    n            = n_subject,
    d            = d,
    fitCOX       = fitCOX,
    prep         = prep,
    # ---- timing / profiling (diagnostic, passed through from RefinedfastEM_LLA) ----
    EMtime       = theta$EMtime,
    EM.time      = theta$EM.time,
    comp.time    = theta$comp.time,
    prof         = theta$prof
  )

  # Store B-spline dataMat for TimeVar models (needed by timevar_eta_plot etc.)
  if (model == "TimeVar") {
    out$dataMat <- list(B = B_list)
    out$model   <- "TimeVar"
    out$Bs      <- Bs
  } else {
    out$model   <- "Fixed"
  }

  out$dataset <- list(
    LMM_dat  = train_dataset$LMM_dat,
    Surv_dat = train_dataset$Surv_dat,
    s = s,
    var_list  = train_dataset$var_list
  )

  class(out) <- "JEL"
  out
}