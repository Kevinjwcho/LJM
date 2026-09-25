############################################################
## pred_surv_compute_LLA.R
## - Survival prediction for LJM fitted objects
## - Uses fitLLAJEL$fitCOX (already stored inside fitLLAJEL)
## - Uses subject-level bLLA (no LME BLUPs)
## - Supports both "Fixed" and "TimeVar" (B-spline eta) models.
##   * Fixed:   length(eta) == 2*nK, eta_b is constant in t.
##   * TimeVar: eta is expanded over B-spline basis per biomarker
##              (nested list stored in fitLLAJEL$dataMat$B). Mixed
##              fits (some biomarkers tv, some fixed via list(NULL))
##              are handled the same way as pred_surv_compute.R's
##              tv branch.
############################################################

pred_surv_compute_LLA <- function(fitLLAJEL,
                                  Cox_dat.te,
                                  b_new,
                                  tau = NULL,
                                  rho = 0,
                                  CI = FALSE,
                                  MC = 100,
                                  alpha = 0.05) {
  # -----------------------------
  # Checks
  # -----------------------------
  if (is.null(fitLLAJEL$fitCOX)) stop("fitLLAJEL must contain $fitCOX.")
  if (is.null(fitLLAJEL$coefficients$hazard)) stop("fitLLAJEL$coefficients$hazard is missing.")
  if (is.null(fitLLAJEL$coefficients$eta)) stop("fitLLAJEL$coefficients$eta is missing.")
  if (is.null(fitLLAJEL$coefficients$phi)) stop("fitLLAJEL$coefficients$phi is missing.")
  if (is.null(fitLLAJEL$prep$s)) stop("fitLLAJEL$prep$s is missing (landmark time).")
  if (!is.data.frame(Cox_dat.te)) stop("Cox_dat.te must be a data.frame (or tibble).")
  if (!("id" %in% names(Cox_dat.te))) stop("Cox_dat.te must contain column: id.")
  if (!all(c("start", "stop", "event") %in% names(Cox_dat.te))) {
    stop("Cox_dat.te must contain columns: start, stop, event (time-dependent Cox format).")
  }

  if (!is.matrix(b_new)) stop("b_new must be a matrix.")
  if (ncol(b_new) %% 2 != 0) stop("b_new must have 2 columns per outcome (b0,b1).")
  if (is.null(colnames(b_new))) stop("b_new must have colnames (e.g., Y.1_b0, Y.1_b1, ...).")
  if (is.null(rownames(b_new))) stop("b_new must have rownames = subject ids (as character).")

  # Detect TimeVar
  is_tv <- identical(fitLLAJEL$model, "TimeVar") && !is.null(fitLLAJEL$dataMat$B)

  # -----------------------------
  # Pull fitted pieces
  # -----------------------------
  fitCOX <- fitLLAJEL$fitCOX

  haz <- fitLLAJEL$coefficients$hazard
  surv_t_lamb <- as.numeric(haz[, 1])
  bashaz      <- as.numeric(haz[, 2])

  eta <- fitLLAJEL$coefficients$eta
  phi <- fitLLAJEL$coefficients$phi
  if (length(phi) == 0) phi <- 0

  # -----------------------------
  # Helper: rho transform
  # -----------------------------
  trans_fun <- function(x) {
    if (rho == 0) x else log(1 + rho * x) / rho
  }

  # -----------------------------
  # Build Cox design for test set (same terms as fitCOX)
  # IMPORTANT: do NOT assume intercept position; follow fitCOX terms.
  # -----------------------------
  TermsSurv <- fitCOX$terms
  formSurv  <- formula(fitCOX)

  mfSurv_te <- model.frame(TermsSurv, Cox_dat.te)

  W_new_full <- model.matrix(formSurv, mfSurv_te)
  if (attr(TermsSurv, "intercept")) {
    W_new <- W_new_full[, -1, drop = FALSE]
  } else {
    W_new <- W_new_full
  }
  n_surv <- nrow(W_new)

  # -----------------------------
  # Align b_new rows to Cox_dat.te (by id)
  # Cox_dat.te can be start/stop rows; we need per-row bLLA
  # -----------------------------
  id_te <- as.character(Cox_dat.te$id)
  idx   <- match(id_te, rownames(b_new))
  if (anyNA(idx)) {
    bad <- unique(id_te[is.na(idx)])
    stop(paste0("b_new is missing these ids: ", paste(bad, collapse = ", ")))
  }
  b_row <- b_new[idx, , drop = FALSE]  # n_surv x (2*nK)

  eta_vec <- as.numeric(eta)
  q <- ncol(b_new)   # = 2*nK

  # -----------------------------
  # Phi vector (shared by both branches)
  # -----------------------------
  phi_vec <- as.numeric(phi)
  if (ncol(W_new) > 0 && length(phi_vec) != ncol(W_new)) {
    stop(sprintf("phi length mismatch. Need %d (ncol(W_new)), got %d.",
                 ncol(W_new), length(phi_vec)))
  }

  base_cumhaz <- cumsum(bashaz)   # used in Fixed branch

  # =============================================================
  # TimeVar branch
  # =============================================================
  if (is_tv) {
    B_list <- fitLLAJEL$dataMat$B
    nK <- q / 2
    if (!is.list(B_list) || length(B_list) != nK) {
      stop(sprintf("fitLLAJEL$dataMat$B must be a list of length nK (=%d) for TimeVar prediction.",
                   nK))
    }

    # Evaluate each basis at surv_t_lamb.
    # - list(NULL)           -> constant 1-column (fixed-eta within tv fit)
    # - bs-class basis obj   -> predict(bs, surv_t_lamb)
    BB <- lapply(seq_along(B_list), function(k) {
      bk <- B_list[[k]]
      if (is.null(unlist(bk))) {
        matrix(1, nrow = length(surv_t_lamb), ncol = 1)
      } else {
        as.matrix(predict(bk, surv_t_lamb))
      }
    })

    eta.inds <- sapply(BB, function(x) 2 * ncol(x))
    if (sum(eta.inds) != length(eta_vec)) {
      stop(sprintf("eta length mismatch in TimeVar path: expected %d (= 2 * sum(df_k)), got %d.",
                   sum(eta.inds), length(eta_vec)))
    }

    # Duplicate each basis column; order matches eta layout produced by
    # InitVal_LLAJEL_tv (b0/b1 interleaved per basis element per biomarker).
    BBc  <- do.call(cbind, BB)                                         # T x sum(d_k)
    BBcd <- do.call(cbind, lapply(seq_len(ncol(BBc)),                  # T x 2*sum(d_k)
                                  function(i) cbind(BBc[, i], BBc[, i])))

    etamat <- matrix(eta_vec, nrow = nrow(BBcd), ncol = ncol(BBcd),
                     byrow = TRUE)
    BBeta  <- BBcd * etamat

    add_odd_even_columns <- function(M) {
      odd  <- seq(1, ncol(M), 2)
      even <- seq(2, ncol(M), 2)
      cbind(rowSums(M[, odd,  drop = FALSE]),
            rowSums(M[, even, drop = FALSE]))
    }

    # etaBBs_mat[t, ] = (eta_b0_k(t), eta_b1_k(t)) stacked over k=1..nK
    etaBBs_mat <- matrix(0, nrow = nrow(BBeta), ncol = 2 * length(eta.inds))
    start <- 1L; end <- 0L
    for (k in seq_along(eta.inds)) {
      end <- end + eta.inds[k]
      etaBB_k  <- BBeta[, start:end, drop = FALSE]
      etaBBs_k <- add_odd_even_columns(etaBB_k)
      etaBBs_mat[, (2*k - 1):(2*k)] <- etaBBs_k
      start <- end + 1L
    }

    # Subject-wise time-varying linear predictor on the hazard grid.
    # eta_b_mat[t, i] = etaBBs_mat[t, ] %*% b_new[i, ]
    eta_b_mat <- etaBBs_mat %*% t(b_row)             # T x n_surv

    if (ncol(W_new) == 0) {
      lp_mat <- eta_b_mat
    } else {
      Wphi <- as.numeric(W_new %*% phi_vec)           # n_surv
      lp_mat <- eta_b_mat +
        matrix(Wphi, nrow = nrow(eta_b_mat), ncol = length(Wphi), byrow = TRUE)
    }

    exp_h_mat <- exp(lp_mat)                          # T x n_surv

    # Per-subject cumulative hazard: cumsum over t of exp_h * bashaz
    ch_mat_tx <- apply(exp_h_mat * bashaz, 2, cumsum) # T x n_surv

    surv_mat <- t(exp(-trans_fun(ch_mat_tx)))         # n_surv x T (align with Fixed branch)

    # -----------------------------
    # Optional CI via Monte Carlo on eta (and phi)
    # -----------------------------
    if (CI) {
      if (is.null(fitLLAJEL$Vcov)) stop("fitLLAJEL$Vcov is missing. Set CI=FALSE or fit with Vcov.")
      Vcov <- fitLLAJEL$Vcov

      eta_ind <- grepl("eta", rownames(Vcov))
      phi_ind <- grepl("phi", rownames(Vcov))
      if (sum(eta_ind) < length(eta_vec)) stop("Could not find eta block in Vcov (TimeVar).")
      if (ncol(W_new) > 0 && sum(phi_ind) < ncol(W_new)) stop("Could not find phi block in Vcov.")

      eta_names <- rownames(Vcov)[eta_ind][seq_along(eta_vec)]
      eta_S     <- Vcov[eta_names, eta_names, drop = FALSE]

      if (ncol(W_new) == 0) {
        phi_S  <- matrix(0, 0, 0)
      } else {
        phi_names <- rownames(Vcov)[phi_ind][seq_len(ncol(W_new))]
        phi_S     <- Vcov[phi_names, phi_names, drop = FALSE]
      }

      eta_MC <- MASS::mvrnorm(MC, mu = eta_vec, Sigma = eta_S)
      if (ncol(W_new) > 0) {
        phi_MC <- MASS::mvrnorm(MC, mu = phi_vec, Sigma = phi_S)
      } else {
        phi_MC <- matrix(0, nrow = MC, ncol = 0)
      }

      surv_array <- array(NA_real_, dim = c(n_surv, length(surv_t_lamb), MC))
      for (mm in seq_len(MC)) {
        eta_mm  <- as.numeric(eta_MC[mm, ])
        etamat_mm <- matrix(eta_mm, nrow = nrow(BBcd), ncol = ncol(BBcd), byrow = TRUE)
        BBeta_mm  <- BBcd * etamat_mm

        etaBBs_mat_mm <- matrix(0, nrow = nrow(BBeta_mm), ncol = 2 * length(eta.inds))
        start <- 1L; end <- 0L
        for (k in seq_along(eta.inds)) {
          end <- end + eta.inds[k]
          etaBB_k  <- BBeta_mm[, start:end, drop = FALSE]
          etaBBs_k <- add_odd_even_columns(etaBB_k)
          etaBBs_mat_mm[, (2*k - 1):(2*k)] <- etaBBs_k
          start <- end + 1L
        }
        eta_b_mm <- etaBBs_mat_mm %*% t(b_row)
        if (ncol(W_new) == 0) {
          lp_mm <- eta_b_mm
        } else {
          Wphi_mm <- as.numeric(W_new %*% as.numeric(phi_MC[mm, ]))
          lp_mm <- eta_b_mm +
            matrix(Wphi_mm, nrow = nrow(eta_b_mm), ncol = length(Wphi_mm), byrow = TRUE)
        }
        ch_mm <- apply(exp(lp_mm) * bashaz, 2, cumsum)
        surv_array[, , mm] <- t(exp(-trans_fun(ch_mm)))
      }

      surv_mean  <- apply(surv_array, c(1, 2), mean)
      surv_lower <- apply(surv_array, c(1, 2), quantile, probs = alpha / 2)
      surv_upper <- apply(surv_array, c(1, 2), quantile, probs = 1 - alpha / 2)

      make_tab <- function(i) {
        data.frame(
          tau = surv_t_lamb,
          pred_surv = surv_mean[i, ],
          lower = surv_lower[i, ],
          upper = surv_upper[i, ]
        )
      }
      surv_tab <- lapply(seq_len(n_surv), make_tab)
    } else {
      make_tab <- function(i) {
        data.frame(
          tau = surv_t_lamb,
          pred_surv = as.numeric(surv_mat[i, ])
        )
      }
      surv_tab <- lapply(seq_len(n_surv), make_tab)
    }

  } else {
    # =============================================================
    # Fixed-eta branch (original behavior)
    # =============================================================
    if (length(eta_vec) != q) {
      stop(sprintf("eta length mismatch. Need %d (=2*nK), got %d.", q, length(eta_vec)))
    }
    eta_b <- as.numeric(b_row %*% eta_vec)             # length n_surv

    if (ncol(W_new) == 0) {
      lp <- eta_b
    } else {
      lp <- as.numeric(W_new %*% phi_vec) + eta_b
    }

    exp_h       <- exp(lp)                              # length n_surv
    cum_haz_mat <- outer(exp_h, base_cumhaz, "*")       # n_surv x T

    surv_mat <- exp(-trans_fun(cum_haz_mat))            # n_surv x T

    if (CI) {
      if (is.null(fitLLAJEL$Vcov)) stop("fitLLAJEL$Vcov is missing. Set CI=FALSE or fit with Vcov/post-process.")
      Vcov <- fitLLAJEL$Vcov

      eta_ind <- grepl("^eta", rownames(Vcov)) | grepl("eta", rownames(Vcov))
      phi_ind <- grepl("^phi", rownames(Vcov)) | grepl("phi", rownames(Vcov))
      if (sum(eta_ind) < q) stop("Could not find eta block in Vcov (name mismatch).")
      if (ncol(W_new) > 0 && sum(phi_ind) < ncol(W_new)) stop("Could not find phi block in Vcov (name mismatch).")

      eta_names <- rownames(Vcov)[eta_ind][seq_len(q)]
      eta_mu    <- eta_vec
      eta_S     <- Vcov[eta_names, eta_names, drop = FALSE]

      if (ncol(W_new) == 0) {
        phi_mu <- numeric(0)
        phi_S  <- matrix(0, 0, 0)
      } else {
        phi_names <- rownames(Vcov)[phi_ind][seq_len(ncol(W_new))]
        phi_mu    <- phi_vec
        phi_S     <- Vcov[phi_names, phi_names, drop = FALSE]
      }

      eta_MC <- MASS::mvrnorm(MC, mu = eta_mu, Sigma = eta_S)
      if (ncol(W_new) > 0) {
        phi_MC <- MASS::mvrnorm(MC, mu = phi_mu, Sigma = phi_S)
      } else {
        phi_MC <- matrix(0, nrow = MC, ncol = 0)
      }

      surv_array <- array(NA_real_, dim = c(n_surv, length(surv_t_lamb), MC))
      for (mm in seq_len(MC)) {
        eta_b_mm <- as.numeric(b_row %*% as.numeric(eta_MC[mm, ]))
        if (ncol(W_new) == 0) {
          lp_mm <- eta_b_mm
        } else {
          lp_mm <- as.numeric(W_new %*% as.numeric(phi_MC[mm, ])) + eta_b_mm
        }
        exp_h_mm   <- exp(lp_mm)
        cum_haz_mm <- outer(exp_h_mm, base_cumhaz, "*")
        surv_array[, , mm] <- exp(-trans_fun(cum_haz_mm))
      }

      surv_mean  <- apply(surv_array, c(1, 2), mean)
      surv_lower <- apply(surv_array, c(1, 2), quantile, probs = alpha / 2)
      surv_upper <- apply(surv_array, c(1, 2), quantile, probs = 1 - alpha / 2)

      make_tab <- function(i) {
        data.frame(
          tau = surv_t_lamb,
          pred_surv = surv_mean[i, ],
          lower = surv_lower[i, ],
          upper = surv_upper[i, ]
        )
      }
      surv_tab <- lapply(seq_len(n_surv), make_tab)
    } else {
      make_tab <- function(i) {
        data.frame(
          tau = surv_t_lamb,
          pred_surv = as.numeric(surv_mat[i, ])
        )
      }
      surv_tab <- lapply(seq_len(n_surv), make_tab)
    }
  }

  # -----------------------------
  # Name by id (row-level ids; if repeated, keep row index)
  # -----------------------------
  nm <- id_te
  if (any(duplicated(nm))) {
    nm <- paste0(nm, "_row", ave(seq_along(nm), nm, FUN = seq_along))
  }
  names(surv_tab) <- nm

  # If tau specified: return one-row-per-id (closest <= tau)
  if (!is.null(tau)) {
    pick1 <- function(df) {
      ii <- which(df$tau <= tau)
      if (length(ii) == 0) return(df[1, , drop = FALSE])
      df[max(ii), , drop = FALSE]
    }
    sub_pred <- lapply(surv_tab, pick1)
    out <- data.frame(id = names(sub_pred), do.call(rbind, sub_pred), row.names = NULL)
    return(out)
  }

  surv_tab
}
