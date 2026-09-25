############################################################
# LJM utilities + main function
# - LLA (kernel-weighted local linear) summaries
# - Cox fit on landmark survival data (optionally include b0/b1)
# - InitVal_LLAJEL() + RefinedfastEM_LLA() pipeline
############################################################

library(dplyr)
library(survival)

# ---------------------------------------------------------
# Kernel weight
# ---------------------------------------------------------
kernel_weight <- function(u, kernel = c("gaussian", "epanechnikov")) {
  kernel <- match.arg(kernel)
  if (kernel == "gaussian") {
    dnorm(u)
  } else {
    w <- 0.75 * (1 - u^2)
    w[abs(u) > 1] <- 0
    w
  }
}

# ---------------------------------------------------------
# Attach LLA summaries to survival data
# Adds columns: <outcome>_b0 and <outcome>_b1
# ---------------------------------------------------------
attach_LLA_to_surv <- function(Surv_dat, bb_list) {
  out <- Surv_dat
  
  y_names <- names(bb_list)
  if (is.null(y_names) || any(y_names == "")) {
    stop("bb_list must be a named list (names(bb_list) = y_vars).")
  }
  if (!("id" %in% names(out))) {
    stop("Surv_dat must contain column: id.")
  }
  
  for (k in seq_along(bb_list)) {
    y  <- y_names[k]
    bb <- bb_list[[k]]
    
    if (!all(c("id", "b0", "b1") %in% names(bb))) {
      stop(paste0("bb_list[['", y, "']] must contain columns: id, b0, b1"))
    }
    
    bb2 <- bb[, c("id", "b0", "b1")]
    names(bb2)[names(bb2) == "b0"] <- paste0(y, "_b0")
    names(bb2)[names(bb2) == "b1"] <- paste0(y, "_b1")
    
    out <- dplyr::left_join(out, bb2, by = "id")
  }
  
  out
}

# ---------------------------------------------------------
# Build bLLA matrices aligned to subject order
# Returns: named list (per outcome) matrix [n_subject x 2]
# ---------------------------------------------------------
make_bLLA_from_list <- function(bb_list, subject_ids, require_ok = FALSE) {
  
  bLLA <- vector("list", length(bb_list))
  names(bLLA) <- names(bb_list)
  
  for (k in seq_along(bb_list)) {
    yname <- names(bb_list)[k]
    df <- bb_list[[k]]
    
    if (require_ok && ("ok" %in% names(df))) {
      df <- df[df$ok %in% TRUE, , drop = FALSE]
    }
    
    mat <- matrix(NA_real_, nrow = length(subject_ids), ncol = 2)
    colnames(mat) <- c(paste0(yname, "_b0"), paste0(yname, "_b1"))
    
    m <- match(subject_ids, df$id)
    good <- !is.na(m)
    
    mat[good, 1] <- df$b0[m[good]]
    mat[good, 2] <- df$b1[m[good]]
    
    bLLA[[k]] <- mat
  }
  
  bLLA
}

# ---------------------------------------------------------
# Build full LLA components (LME-like structure)
# ---------------------------------------------------------
make_LLA_full <- function(LMM_dat, y_vars, s, h,
                          ker = "gaussian",
                          drop_zero_weight = TRUE,
                          drop_na_y = TRUE,
                          subject_ids = NULL,
                          var_list = list(id = "id", time = "time")) {
  
  LLA_summary <- vector("list", length(y_vars))
  Y_list <- vector("list", length(y_vars))
  Z_list <- vector("list", length(y_vars))
  W_list <- vector("list", length(y_vars))
  ID_list <- vector("list", length(y_vars))
  ni_list <- vector("list", length(y_vars))
  uniqueID_list <- vector("list", length(y_vars))
  
  names(LLA_summary) <- names(Y_list) <- names(Z_list) <-
    names(W_list) <- names(ID_list) <- names(ni_list) <-
    names(uniqueID_list) <- y_vars
  
  for (k in seq_along(y_vars)) {
    
    yname <- y_vars[k]
    
    id   <- LMM_dat[[var_list[["id"]]]]
    time <- LMM_dat[[var_list[["time"]]]]
    y    <- LMM_dat[[yname]]
    
    u <- (time - s) / h
    w <- kernel_weight(u, kernel = ker)
    
    ok <- rep(TRUE, length(y))
    if (drop_na_y) ok <- ok & !is.na(y)
    if (drop_zero_weight) ok <- ok & (w > 0)
    
    id   <- id[ok]
    time <- time[ok]
    y    <- y[ok]
    w    <- w[ok]

    # ------------------------------------------------------------------
    # Keep every survival subject present in THIS outcome's stack. A
    # bounded kernel (epanechnikov / uniform) at a narrow bandwidth can
    # zero-weight ALL of a subject's observations for this outcome, so
    # the drop_zero_weight step above removes the subject entirely. That
    # makes this outcome's per-subject lists (Z.st/Y.st/Wker.st/ni)
    # shorter than n_subject, and the EM's positional indexing
    # Z.st[[k]][[i]] (RefinedfastEM_LLA.R) then runs out of bounds.
    # Re-insert any missing subject with a single placeholder row
    # (time = s -> x = 0 -> Z = (1, 0), y = 0).
    #
    # Weight: use a TINY positive w_ph (not exactly 0). A hard 0 makes the
    # subject's kernel-weighted precision diagonal exactly 0, and the C++
    # E-step rejects it ("inv_omega diagonal must be positive"). A tiny
    # positive weight keeps that diagonal strictly positive while adding a
    # negligible amount (~w_ph) to the weighted LS / residual variance /
    # Bsigma; the subject's local LLA is still NA -> imputed with cLLA
    # (the same "zero deviation from the population mean" used elsewhere).
    # ------------------------------------------------------------------
    if (!is.null(subject_ids)) {
      miss <- setdiff(subject_ids, unique(id))
      if (length(miss) > 0) {
        w_ph <- 1e-6
        id   <- c(id,   miss)
        time <- c(time, rep(s,   length(miss)))
        y    <- c(y,    rep(0,    length(miss)))
        w    <- c(w,    rep(w_ph, length(miss)))
      }
    }

    ord <- order(id, time)
    id   <- id[ord]
    time <- time[ord]
    y    <- y[ord]
    w    <- w[ord]
    
    x <- time - s
    Z <- cbind(1, x)
    colnames(Z) <- c("Intercept", "time_shift")
    
    uniqueID <- !duplicated(id)
    idx <- which(uniqueID)
    ni  <- diff(c(idx, length(id) + 1))
    ID  <- rep(seq_along(ni), times = ni)
    
    Y_list[[k]] <- as.vector(y)
    Z_list[[k]] <- Z
    W_list[[k]] <- as.vector(w)
    ID_list[[k]] <- ID
    ni_list[[k]] <- ni
    uniqueID_list[[k]] <- uniqueID
    
    df_tmp <- data.frame(id = id, time = time, y = y, w = w)
    
    summary_df <- df_tmp %>%
      dplyr::group_by(id) %>%
      dplyr::group_modify(~{
        d <- .x
        neff <- sum(d$w, na.rm = TRUE)
        
        if (nrow(d) < 2 || neff < 1e-8) {
          return(data.frame(b0 = NA_real_, b1 = NA_real_, n_eff = neff, ok = FALSE))
        }
        
        x <- d$time - s
        yv <- d$y
        wv <- d$w
        
        X <- cbind(1, x)
        XtW <- t(X) * wv
        XtWX <- XtW %*% X
        XtWy <- XtW %*% yv
        
        if (det(XtWX) < 1e-12) {
          return(data.frame(b0 = NA_real_, b1 = NA_real_, n_eff = neff, ok = FALSE))
        }
        
        coef <- solve(XtWX, XtWy)
        
        data.frame(
          b0 = as.numeric(coef[1]),
          b1 = as.numeric(coef[2]),
          n_eff = neff,
          ok = TRUE
        )
      }) %>%
      dplyr::ungroup()
    
    LLA_summary[[k]] <- summary_df
  }
  
  list(
    LLA_summary = LLA_summary,
    Y = Y_list,
    Z = Z_list,
    W = W_list,
    ID = ID_list,
    ni = ni_list,
    uniqueID = uniqueID_list
  )
}

# ---------------------------------------------------------
# Bandwidth that places EVERY longitudinal observation inside the kernel
# support, so the initial values are built from the FULL longitudinal record
# with the kernel's own shape (largest weight at the landmark, tapering to the
# edges) rather than from a localized sub-window. Compact-support kernels
# (Epanechnikov is zero at |u| = 1) need h strictly greater than the observed
# range, hence the (1 + eps) factor; the Gaussian kernel has infinite support
# so the range alone already retains every observation. The kernel itself is
# always `ker`, i.e. the SAME kernel option the EM fit uses -- only the
# bandwidth differs, and it is not user-tunable.
# ---------------------------------------------------------
# FULL_RANGE_PAD: how far past the observed range the initial-value bandwidth
# is stretched. Only matters for compact-support kernels: with pad = 1 the
# Epanechnikov weight at the far edge is exactly 0, and a pad barely above 1
# makes it ~0, which `varFixed(~ 1/w)` turns into a near-infinite residual
# variance for that observation. 1.05 keeps every observation at a usable
# weight (edge/centre ratio about 11) while still covering the whole record.
# The Gaussian kernel is insensitive to this (infinite support).
FULL_RANGE_PAD <- 1.05

full_range_h <- function(times, s, pad = FULL_RANGE_PAD) {
  R <- max(abs(as.numeric(times) - s), na.rm = TRUE)
  if (!is.finite(R) || R <= 0) return(1)
  R * pad
}

# ---------------------------------------------------------
# Pooled initialization via the kernel-weighted local linear mixed model.
# For each biomarker, fit  y ~ 1 + (t - s)  with random intercept & slope and
# residual variance sigma^2 / K(t - s) (varFixed on 1/weight), i.e. the
# longitudinal sub-model of Section 2.1 over the FULL longitudinal record
# (bandwidth full_range_h, kernel `ker`, no observation excluded). This is the pooled
# ("kernel LMME") counterpart of the per-subject WLS start: pooling via
# b_i ~ N(c, Sigma_b) stabilizes sparse subjects, and the model's own REML
# estimates give consistent, non-shrunk starting values.
# Returns (all aligned to subject_ids / stacked by (b0,b1) per marker):
#   bLLA    : list of K (n_subject x 2) BLUP matrices (value & slope at s)
#   cLLA    : list of K length-2 fixed effects  = c
#   Bsigma  : 2K x 2K block-diagonal Sigma_b(s) from getVarCov (cross-marker 0)
#   Ysigma2 : length-K residual variances sigma^2
# NULL-safe: caller falls back to per-subject WLS init if this errors.
# ---------------------------------------------------------
make_pooled_init <- function(LMM_dat, y_vars, s, ker, subject_ids,
                             var_list = list(id = "id", time = "time"),
                             h_init = NULL) {
  # h_init: bandwidth of the pooled kernel LMM (JEL 2.4: the EM bandwidth h by
  # default; the caller falls back to full_range_h() if this fit fails).
  if (is.null(h_init)) h_init <- full_range_h(LMM_dat[[var_list[["time"]]]], s)
  idv <- var_list[["id"]]; tv <- var_list[["time"]]
  sid <- as.character(subject_ids); nS <- length(sid); K <- length(y_vars)
  bLLA <- vector("list", K); names(bLLA) <- y_vars
  cLLA <- vector("list", K)
  Ysigma2 <- numeric(K)
  Bsigma  <- matrix(0, 2L * K, 2L * K)
  for (k in seq_len(K)) {
    tc <- as.numeric(LMM_dat[[tv]]) - s
    w  <- kernel_weight(tc / h_init, kernel = ker)  # same kernel as the EM fit, at h_init
    y  <- as.numeric(LMM_dat[[y_vars[k]]])
    id <- as.character(LMM_dat[[idv]])
    ok <- !is.na(y) & is.finite(w) & (w > 0)
    d  <- data.frame(id = droplevels(factor(id[ok], levels = sid)),
                     tc = tc[ok], y = y[ok], iw = 1 / w[ok],
                     stringsAsFactors = FALSE)
    ctl <- nlme::lmeControl(opt = "optim", returnObject = TRUE,
                            msMaxIter = 200, msMaxEval = 400)
    # Diagonal random effects (independent intercept & slope) are far more
    # stable under kernel down-weighting than an unstructured 2x2; the
    # intercept-slope covariance is left 0 in the initial value and filled in
    # by the EM. Fall back to random-intercept-only, then error (-> WLS init).
    fit <- tryCatch(
      nlme::lme(y ~ tc, random = list(id = nlme::pdDiag(~ tc)), data = d,
                weights = nlme::varFixed(~ iw), method = "REML", control = ctl),
      error = function(e)
        nlme::lme(y ~ tc, random = ~ 1 | id, data = d,
                  weights = nlme::varFixed(~ iw), method = "REML", control = ctl))
    cf <- nlme::fixef(fit)                   # (Intercept), tc = value & slope at s
    co <- stats::coef(fit)                   # per-subject BLUPs
    Graw <- as.matrix(nlme::getVarCov(fit))  # 1x1 (int only) or 2x2
    G <- matrix(0, 2, 2)
    if (all(dim(Graw) == c(2,2))) { G <- Graw } else { G[1,1] <- Graw[1,1]; G[2,2] <- var(co[["tc"]], na.rm=TRUE) }
    M  <- matrix(rep(cf, each = nS), nrow = nS)   # default rows = population mean
    idx <- match(rownames(co), sid); okk <- !is.na(idx)
    M[idx[okk], 1] <- co[["(Intercept)"]][okk]
    if ("tc" %in% names(co)) M[idx[okk], 2] <- co[["tc"]][okk]  # else keep population slope
    colnames(M) <- c("b0", "b1")
    bLLA[[k]] <- M
    cLLA[[k]] <- as.numeric(cf)
    Bsigma[(2L * k - 1L):(2L * k), (2L * k - 1L):(2L * k)] <- G
    Ysigma2[k] <- fit$sigma^2
  }
  list(bLLA = bLLA, cLLA = cLLA, Bsigma = Bsigma, Ysigma2 = Ysigma2)
}

# ---------------------------------------------------------
# Landmark preparation (returns prep list)
# ---------------------------------------------------------
prep_LLA_landmark <- function(LMM_dat, Surv_dat, y_vars, s, h, ker = "gaussian",
                              var_list = list(id = "id", time = "time"),
                              init_pooled = TRUE) {

  # h      : FITTING bandwidth -> EM kernel weights (Y/Z/W_ker) and residual var.
  # Initial values (bLLA/cLLA/Bsigma, Cox init) are ALWAYS built from the full
  # longitudinal record at full_range_h(), using the same kernel `ker` as the
  # EM fit. The kernel shape is preserved, no observation is dropped, and there
  # is no separate init-bandwidth argument; the EM then fits at the user's h.
  # NOTE: the >= min_window_obs subject filter lives in LJM_dat() so that the
  #       whole train_dataset (LMM_dat + Surv_dat + all downstream indexing) is
  #       filtered consistently; do not filter here.

  subject_ids <- Surv_dat$id[!duplicated(Surv_dat$id)]

  # EM data + kernel weights at the FITTING bandwidth h
  LLA_full <- make_LLA_full(
    LMM_dat = LMM_dat,
    y_vars  = y_vars,
    s      = s,
    h       = h,
    ker     = ker,
    subject_ids = subject_ids,
    var_list = var_list
  )

  # Initial-value summary over the full record (same kernel, nothing dropped)
  h_full <- full_range_h(LMM_dat[[var_list[["time"]]]], s)
  LLA_init <- make_LLA_full(
    LMM_dat = LMM_dat,
    y_vars  = y_vars,
    s      = s,
    h       = h_full,
    ker     = ker,
    subject_ids = subject_ids,
    var_list = var_list
  )

  LLA_list_s <- LLA_init$LLA_summary
  
  surv_s2 <- attach_LLA_to_surv(
    Surv_dat = Surv_dat,
    bb_list  = LLA_list_s
  )
  
  bLLA <- make_bLLA_from_list(
    bb_list     = LLA_list_s,
    subject_ids = subject_ids,
    require_ok  = FALSE
  )

  # Pooled (kernel LMME) initialization: fit the longitudinal sub-model directly
  # and take c, Sigma_b(s), sigma^2 and the BLUPs from its REML estimates. When
  # this succeeds we use those values as-is and skip the moment-based summaries
  # below (which underestimate Sigma_b from shrunk BLUPs).
  pooled_init <- NULL; h_init_used <- NA_real_; init_mode <- "wls_full_range"
  if (isTRUE(init_pooled)) {
    # JEL 2.4 rule: the pooled kernel LMM is fitted at the SAME bandwidth h as the
    # EM (pooling across subjects stabilises sparse subjects, so no wider window
    # is needed). If that fit fails, or returns a degenerate variance component,
    # fall back to the pooled fit over the full record (the 2.3 rule); only if
    # that fails too do we use the per-subject WLS start.
    .degenerate <- function(pi) !all(is.finite(pi$Ysigma2)) || any(pi$Ysigma2 <= 0) ||
                                !all(is.finite(diag(pi$Bsigma))) || any(diag(pi$Bsigma) <= 1e-10)
    for (cand in list(list(h = h, mode = "pooled_h"), list(h = h_full, mode = "pooled_full_range"))) {
      pi <- tryCatch(make_pooled_init(LMM_dat, y_vars, s, ker, subject_ids, var_list, h_init = cand$h),
                     error = function(e) { message("prep_LLA_landmark: pooled init at h=",
                       signif(cand$h, 3), " failed (", conditionMessage(e), ")"); NULL })
      if (!is.null(pi) && .degenerate(pi)) {
        message("prep_LLA_landmark: pooled init at h=", signif(cand$h, 3), " degenerate (zero variance component); trying next")
        pi <- NULL
      }
      if (!is.null(pi)) { pooled_init <- pi; h_init_used <- cand$h; init_mode <- cand$mode; break }
    }
    if (is.null(pooled_init)) message("prep_LLA_landmark: all pooled inits failed; using per-subject WLS init.")
    if (!is.null(pooled_init)) bLLA <- pooled_init$bLLA
  }

  # ------------------------------------------------------------------
  # NA-imputation for subjects whose kernel-weighted local LLA failed
  # (e.g., effective sample size near zero, X'WX singular).
  #   - JEL <= 1.2 sidestepped this by pre-windowing data via past_int.
  #   - JEL 2.0 keeps all longitudinal obs and lets the kernel down-
  #     weight far ones, which can leave a subject with no usable obs.
  # Without imputation, NA rows in bLLA propagate to Bsigma (-> all
  # NaN), to D (init random-effects covariance), and the C++ EM fails
  # at inv_sympd(D).  Replacing NA rows with the population center cLLA
  # is equivalent to "this subject has zero deviation from the mean":
  # they contribute 0 to Bsigma, 0 to the Cox association term, and
  # don't bias the rest of the fit.
  # ------------------------------------------------------------------
  na_subjects <- lapply(bLLA, function(mat) which(!complete.cases(mat)))
  cLLA_valid  <- lapply(bLLA, function(mat) colMeans(mat, na.rm = TRUE))

  for (k in seq_along(bLLA)) {
    bad <- na_subjects[[k]]
    if (length(bad) > 0) {
      if (any(!is.finite(cLLA_valid[[k]]))) {
        stop(sprintf(
          "prep_LLA_landmark: outcome '%s' has < 2 subjects with valid local LLA estimates -- the kernel window is too narrow / biomarker is too sparse at this landmark.",
          y_vars[k]
        ))
      }
      bLLA[[k]][bad, ] <- matrix(cLLA_valid[[k]],
                                 nrow  = length(bad),
                                 ncol  = ncol(bLLA[[k]]),
                                 byrow = TRUE)
    }
  }
  n_imputed <- sum(vapply(na_subjects, length, integer(1)))
  if (n_imputed > 0) {
    message(sprintf(
      "prep_LLA_landmark: imputed %d subject-outcome bLLA rows with cLLA (kernel had no usable obs).",
      n_imputed
    ))
  }

  cLLA <- cLLA_valid

  nLLA <- length(y_vars)

  # kernel-weighted residual variance per outcome
  Ysigma2 <- sapply(seq_len(nLLA), function(k){
    Yk <- LLA_full$Y[[k]]
    Wk <- LLA_full$W[[k]]
    Zk <- LLA_full$Z[[k]]
    bk <- bLLA[[k]]
    IDk <- LLA_full$ID[[k]]

    resid <- Yk - rowSums(Zk * bk[IDk, , drop = FALSE])
    sum(Wk * resid^2) / sum(Wk)
  })

  # covariance of (b - c) across outcomes, stacked by (b0,b1).
  # Imputed rows contribute exactly 0 to the cross-product, and we use
  # n_valid - 1 as denominator so the variance estimate isn't shrunk.
  b_minus_c_list <- mapply(function(b, c) sweep(b, 2, c, "-"), bLLA, cLLA, SIMPLIFY = FALSE)
  b_minus_c <- do.call(cbind, b_minus_c_list)
  n_valid_per_outcome <- vapply(na_subjects,
                                function(idx) length(subject_ids) - length(idx),
                                integer(1))
  n_valid_min <- max(min(n_valid_per_outcome), 1L)
  Bsigma <- t(b_minus_c) %*% b_minus_c / max(n_valid_min - 1L, 1L)

  # Prefer the kernel-LMME estimates for the initial c, Sigma_b and sigma^2.
  if (!is.null(pooled_init)) {
    cLLA    <- pooled_init$cLLA
    Bsigma  <- pooled_init$Bsigma
    Ysigma2 <- pooled_init$Ysigma2
  }
  
  list(
    y_vars = y_vars,
    s = s, h = h, ker = ker,
    h_full = h_full,          # full-range bandwidth (fallback for the initial values)
    h_init_used = h_init_used, # bandwidth actually used for the pooled initial values
    init_mode = init_mode,     # pooled_h | pooled_full_range | wls_full_range
    subject_ids = subject_ids,
    LLA_list_s = LLA_list_s,
    surv_s2 = surv_s2,
    bLLA = bLLA,
    cLLA = cLLA,
    Y = LLA_full$Y,
    Z = LLA_full$Z,
    W_ker = LLA_full$W,
    Ysigma2 = Ysigma2,
    Bsigma = Bsigma,
    ID = LLA_full$ID,
    ni = LLA_full$ni,
    uniqueID = LLA_full$uniqueID
  )
}

# ---------------------------------------------------------
# Cox formula builder
# include_LLA = TRUE -> add Y.k_b0 and Y.k_b1 into Cox
# ---------------------------------------------------------
make_cox_formula_from_LLA <- function(LLA_list, base_terms = NULL) {
  
  rhs_terms <- base_terms
  
  # normalise base_terms given as NULL, character(0) or ""
  if (is.null(rhs_terms)) rhs_terms <- character(0)
  rhs_terms <- rhs_terms[!is.na(rhs_terms) & nzchar(rhs_terms)]
  
  rhs <- if (length(rhs_terms) == 0) {
    "1"                 # no baseline covariates
  } else {
    paste(rhs_terms, collapse = " + ")
  }
  
  stats::as.formula(paste0("survival::Surv(start, stop, event) ~ ", rhs))
}