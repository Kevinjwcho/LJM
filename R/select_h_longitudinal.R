#' Longitudinal LSCV bandwidth selection for LJM
#'
#' Selects the kernel bandwidth \code{h} for the Local Linear Approximation
#' (LLA) step of \code{LJM} by leave-one-observation-out cross-validation on
#' the longitudinal data, evaluated via the PRESS leverage identity (no
#' obs-by-obs refit required).
#'
#' The criterion is GLOBAL: every longitudinal observation enters the LOO-CV
#' (there is no window around the landmark), so candidate bandwidths are compared
#' on the full observation set.  Population covariance \code{D} and noise variance
#' \code{var.e} are taken from a single baseline \code{LJM} fit at
#' \code{baseline_h} (default \code{median(h_grid)}); the centering \code{c}
#' (\code{cLLA}) is refreshed at every candidate \code{h}.  For each subject
#' \eqn{i} and candidate \code{h}, the BLUP precision is
#' \deqn{\Lambda_i = D^{-1} + Z_i^\top R_i^{-1} Z_i,
#'       \quad R_i = \mathrm{diag}(\sigma_k^2 / w_{ij}),}
#' the posterior mean is
#' \deqn{\mu_i = c + \Lambda_i^{-1} Z_i^\top R_i^{-1} (Y_i - Z_i c),}
#' and the leave-one-out residual at obs \eqn{(j, k)} follows from
#' \deqn{e_{LOO,jk} = \frac{Y_{ijk} - z_{jk}^\top \mu_i}{1 - h_{jk}},
#'       \quad h_{jk} = \frac{w_{ij}}{\sigma_k^2}\, z_{jk}^\top \Lambda_i^{-1} z_{jk}.}
#' For obs outside the kernel support at this \code{h} (i.e. \eqn{w_{ij} = 0}),
#' \eqn{h_{jk} = 0} and the LOO residual reduces to the regular BLUP residual.
#'
#' Per-process selection minimises \eqn{\mathrm{MSE}_k(h)} separately; the shared
#' criterion is \eqn{\sum_k \mathrm{MSE}_k(h) / \mathrm{var.e}_k}, normalised by
#' baseline-fit noise scale so processes with different magnitudes are
#' comparable.
#'
#' @param LMM_dat Long-format longitudinal data frame containing the columns
#'   named in \code{var_list} plus the outcome columns listed in \code{y_vars}.
#' @param Surv_dat Survival data frame (subject-level) containing the subject id
#'   column named in \code{var_list}.  Used to enumerate the subject roster.
#' @param y_vars Character vector of longitudinal outcome column names.
#' @param s Landmark time (scalar).
#' @param h_grid Numeric vector of candidate bandwidths (length >= 2, all > 0).
#' @param ker Kernel: one of \code{"epanechnikov"} (default) or \code{"gaussian"}.
#' @param var_list Named list with at least \code{id} and \code{time} pointing to
#'   the corresponding column names in \code{LMM_dat} / \code{Surv_dat}.
#' @param baseline_h Bandwidth used for the one-shot baseline EM fit
#'   (\code{LJM}).  Defaults to \code{median(h_grid)}.
#' @param fit Optional pre-fitted \code{"LJM"} object (the output of
#'   \code{\link{LJM}}) at \code{baseline_h}.  If supplied, the baseline EM is
#'   skipped and \code{D}, \code{var.e} are taken from this object.
#' @param train_dataset Training-dataset list of the form expected by
#'   \code{\link{LJM}}.  Required when \code{fit} is \code{NULL}.
#' @param base_terms,Bs,gh.nodes,max.iter,tol Passed through to
#'   \code{\link{LJM}} when running the baseline EM.
#' @param verbose Logical; if \code{TRUE}, prints per-h CV values during the
#'   sweep.
#' @param ... Ignored; absorbs deprecated arguments (e.g. the former \code{W}
#'   criterion window) so existing callers do not error.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{h_grid}}{The input candidate grid.}
#'   \item{\code{cv_per_process}}{Matrix \eqn{|h\_grid| \times nK} of per-process
#'     LOO MSEs.}
#'   \item{\code{cv_shared}}{Numeric vector of scale-normalised shared
#'     criterion values.}
#'   \item{\code{h_per_process}}{Length-\eqn{nK} vector of per-process
#'     minimisers.}
#'   \item{\code{h_shared}}{Single minimiser of \code{cv_shared}.}
#'   \item{\code{n_in_W_per_proc}}{Per-process obs counts inside the criterion
#'     window (identical across rows by construction).}
#'   \item{\code{W}}{Always \code{Inf} (criterion is global; window removed).}
#'   \item{\code{ker}}{Kernel used.}
#'   \item{\code{baseline_h}}{Bandwidth at which the baseline EM was run.}
#'   \item{\code{baseline_fit}}{The baseline \code{"LJM"} fit object.}
#' }
#'
#' @seealso \code{\link{LJM}}
#'
#' @export
select_h_longitudinal <- function(
    LMM_dat,
    Surv_dat,
    y_vars,
    s,
    h_grid,
    ker            = c("epanechnikov", "gaussian"),
    var_list       = list(id = "id", time = "time"),
    baseline_h     = NULL,
    fit            = NULL,
    train_dataset  = NULL,
    base_terms     = character(0),
    Bs             = NULL,
    gh.nodes       = 7,
    max.iter       = 1000,
    tol            = 0.001,
    verbose        = FALSE,
    ...
) {
  ker <- match.arg(ker)
  stopifnot(length(h_grid) >= 2, all(h_grid > 0))

  # ---- 1. Criterion window (REMOVED) ----------------------------------------
  # The LOO-CV / PRESS criterion is now GLOBAL: it uses every longitudinal
  # observation, not a window |t - s| <= W around the landmark. The `W`
  # argument is retained only as a no-op for backward compatibility with
  # existing callers that still pass `W = ...`; its value is ignored.

  # ---- 2. Baseline EM fit (one-time) ----------------------------------------
  if (is.null(fit)) {
    if (is.null(train_dataset))
      stop("`select_h_longitudinal`: provide either `fit` (pre-fitted at ",
           "baseline_h) or `train_dataset` to run a one-shot baseline EM.")
    if (is.null(baseline_h)) baseline_h <- stats::median(h_grid)
    if (verbose) message(sprintf("[select_h] baseline EM at h = %.3f ...",
                                 baseline_h))
    fit <- LJM(
      train_dataset = train_dataset,
      y_vars        = y_vars,
      s             = s,
      h             = baseline_h,
      base_terms    = base_terms,
      ker           = ker,
      Bs            = Bs,
      gh.nodes      = gh.nodes,
      max.iter      = max.iter,
      tol           = tol,
      Vcov          = FALSE,
      verbose       = FALSE
    )
  } else if (is.null(baseline_h)) {
    baseline_h <- tryCatch(fit$h, error = function(e) NA_real_)
    if (is.null(baseline_h)) baseline_h <- NA_real_
  }

  D_fit <- fit$coefficients$D
  vare  <- fit$coefficients$var.e
  nK    <- length(y_vars)
  stopifnot(nrow(D_fit) == 2L * nK, length(vare) == nK)

  D_inv <- tryCatch(
    solve(D_fit),
    error = function(e) chol2inv(chol(D_fit + diag(1e-8, 2L * nK)))
  )

  # ---- 3. Pre-compute criterion obs set F (ALL observations; global CV) ------
  id_col   <- var_list$id
  time_col <- var_list$time
  subject_ids <- sort(unique(Surv_dat[[id_col]]))
  nsub <- length(subject_ids)

  F_list <- vector("list", nK)
  for (k in seq_len(nK)) {
    yname <- y_vars[k]
    ldat  <- data.frame(
      id      = LMM_dat[[id_col]],
      t_shift = LMM_dat[[time_col]] - s,
      y       = LMM_dat[[yname]]
    )
    ldat <- ldat[!is.na(ldat$y), , drop = FALSE]
    ldat <- ldat[order(ldat$id, ldat$t_shift), , drop = FALSE]
    F_list[[k]] <- split(ldat[, c("t_shift", "y"), drop = FALSE], ldat$id)
  }

  # ---- 4. Per-h LSCV via PRESS leverage -------------------------------------
  n_h <- length(h_grid)
  cv_per_process <- matrix(NA_real_, n_h, nK)
  n_in_W         <- matrix(0L,       n_h, nK)

  for (i_h in seq_along(h_grid)) {
    h_cand <- h_grid[i_h]

    prep <- prep_LLA_landmark(
      LMM_dat  = LMM_dat,
      Surv_dat = Surv_dat,
      y_vars   = y_vars,
      s       = s,
      h        = h_cand,
      ker      = ker,
      var_list = var_list
    )
    cvec <- as.numeric(unlist(prep$cLLA))

    sse_per_proc <- numeric(nK)
    n_per_proc   <- integer(nK)

    for (i in seq_len(nsub)) {
      sid_str <- as.character(subject_ids[i])

      Zi_blocks <- vector("list", nK)
      Yi_blocks <- vector("list", nK)
      wi_blocks <- vector("list", nK)
      m_ik      <- integer(nK)

      for (k in seq_len(nK)) {
        subdat <- F_list[[k]][[sid_str]]
        if (is.null(subdat) || nrow(subdat) == 0L) {
          Zi_blocks[[k]] <- matrix(numeric(0), 0L, 2L)
          Yi_blocks[[k]] <- numeric(0)
          wi_blocks[[k]] <- numeric(0)
          m_ik[k]        <- 0L
        } else {
          t_shift <- subdat$t_shift
          Zi_blocks[[k]] <- cbind(1, t_shift)
          Yi_blocks[[k]] <- subdat$y
          wi_blocks[[k]] <- kernel_weight(t_shift / h_cand, kernel = ker)
          m_ik[k]        <- length(t_shift)
        }
      }

      if (sum(m_ik) == 0L) next

      Zi <- as.matrix(Matrix::bdiag(Zi_blocks))
      Yi <- unlist(Yi_blocks)
      wi <- unlist(wi_blocks)
      proc_idx <- rep(seq_len(nK), times = m_ik)
      sig_rep  <- vare[proc_idx]

      Wvec   <- wi / sig_rep
      ZtWZ   <- crossprod(Zi, Zi * Wvec)
      Lambda <- D_inv + ZtWZ
      LamInv <- tryCatch(
        solve(Lambda),
        error = function(e)
          chol2inv(chol(Lambda + diag(1e-8, nrow(Lambda))))
      )

      resid_c <- Yi - as.numeric(Zi %*% cvec)
      score   <- crossprod(Zi, Wvec * resid_c)
      mu      <- cvec + as.numeric(LamInv %*% score)

      yhat       <- as.numeric(Zi %*% mu)
      resid_full <- Yi - yhat
      H_diag     <- Wvec * rowSums((Zi %*% LamInv) * Zi)
      resid_loo  <- resid_full / pmax(1 - H_diag, 1e-8)

      for (k in seq_len(nK)) {
        idx_k <- which(proc_idx == k)
        if (length(idx_k) > 0L) {
          sse_per_proc[k] <- sse_per_proc[k] + sum(resid_loo[idx_k]^2)
          n_per_proc[k]   <- n_per_proc[k]   + length(idx_k)
        }
      }
    }

    n_per_safe <- pmax(n_per_proc, 1L)
    cv_per_process[i_h, ] <- sse_per_proc / n_per_safe
    n_in_W[i_h, ]         <- n_per_proc

    if (verbose) {
      cat(sprintf(
        "[select_h] h = %.3f   CV_k = %s   (n_in_W = %s)\n",
        h_cand,
        paste(sprintf("%.4g", cv_per_process[i_h, ]), collapse = ", "),
        paste(n_per_proc, collapse = ", ")
      ))
    }
  }

  # ---- 5. Aggregate criteria ------------------------------------------------
  cv_shared      <- as.numeric(cv_per_process %*% (1 / vare))
  h_shared       <- h_grid[which.min(cv_shared)]
  h_per_process  <- h_grid[apply(cv_per_process, 2, which.min)]

  list(
    h_grid          = h_grid,
    cv_per_process  = cv_per_process,
    cv_shared       = cv_shared,
    h_per_process   = h_per_process,
    h_shared        = h_shared,
    n_in_W_per_proc = n_in_W,
    W               = Inf,   # global criterion (window removed); W arg ignored
    ker             = ker,
    baseline_h      = baseline_h,
    baseline_fit    = fit
  )
}
