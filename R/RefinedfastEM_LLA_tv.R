# ============================================================
# Refined Fast EM for Kernel-weighted LJM (Time-Varying)
# ------------------------------------------------------------
# Combines:
#   - RefinedfastEM_LLA (kernel-weighted LLA, no X/beta, centered prior)
#   - RefinedfastEM_tv  (B-spline time-varying association eta)
#
# Longitudinal part: Kernel-weighted LLA (unchanged from RefinedfastEM_LLA)
# Survival part: Time-varying eta via B-spline basis (from RefinedfastEM_tv)
#
# Key differences from RefinedfastEM_LLA:
#   - eta is expanded to n_eta = sum(ncol(B_list[[k]]) * 2) coefficients
#   - B_list, eta.inds, BBi, etaBBs used throughout
#   - ll_lla_tv / gradll_lla_tv / sdll_lla_tv replace ll_lla / gradll_lla / sdll_lla
#   - calc_mu_surv_noX_tv / calc_Sigma2_surv_noX_tv replace fixed versions
#   - Esurv_exp_t / Setaphi_lla_tv / Hetaphi_lla_tv replace fixed versions
#   - update_lambda_tv replaces update_lambda
# ============================================================

RefinedfastEM_LLA_tv <- function(
    data,
    n_LLA, ID_surv, ni_surv,
    Z.st, Y.st, Wker.st, Wtime,
    Ysigma2, Bsigma,
    bLLA, cLLA, n_subject,
    surv.init, ni, ncw,
    B_list, eta_n,
    gh.nodes = 3, collect.hist = TRUE, max.iter = 200,
    tol = 0.01, diff.type = "abs.rel",
    post.process = FALSE, verbose = FALSE,
    update_c = TRUE
){

  # ----------------------------------------------------------
  # Basic setup
  # ----------------------------------------------------------
  start.time <- proc.time()[3]

  nK <- n_LLA
  q  <- 2 * nK
  uids <- ID_surv
  n <- length(uids)

  diff   <- 100
  b.diff <- 100
  iter   <- 0

  # ----------------------------------------------------------
  # Construct subject-wise stacked longitudinal objects
  #   (identical to RefinedfastEM_LLA)
  # ----------------------------------------------------------
  Z <- lapply(seq_len(length(ni_surv)), function(i) {
    blocks <- sapply(Z.st, function(sublist) sublist[[i]], simplify = FALSE)
    as.matrix(Matrix::bdiag(blocks))
  })

  Y <- lapply(seq_len(length(ni_surv)), function(i) {
    sub <- sapply(Y.st, function(sublist) sublist[[i]], simplify = FALSE)
    matrix(unlist(sub), ncol = 1)
  })

  mi <- lapply(seq_len(length(ni_surv)), function(i) {
    sapply(ni, function(sublist) sublist[[i]], simplify = TRUE)
  })

  W_kernel <- lapply(seq_len(length(ni_surv)), function(i) {
    sub <- sapply(Wker.st, function(sublist) sublist[[i]], simplify = FALSE)
    as.numeric(unlist(sub))
  })

  Zk <- lapply(seq_len(length(ni_surv)), function(i) {
    sapply(Z.st, function(sublist) sublist[[i]], simplify = FALSE)
  })

  # univariate-safe (nK == 1 with m_i == 1 subjects):
  # never let sapply(simplify=TRUE) collapse to a scalar/vector;
  # always force an m_i x nK matrix so ncol() / [, k] work.
  Yk <- lapply(seq_len(length(ni_surv)), function(i) {
    sub <- lapply(Y.st, function(sublist) as.numeric(sublist[[i]]))
    do.call(cbind, sub)
  })

  Wkerk <- lapply(seq_len(length(ni_surv)), function(i) {
    sub <- lapply(Wker.st, function(sublist) as.numeric(sublist[[i]]))
    do.call(cbind, sub)
  })

  w_longK <- lapply(Wkerk, function(M) {
    M <- as.matrix(M)
    lapply(seq_len(ncol(M)), function(k) as.numeric(M[, k]))
  })

  K <- lapply(seq_len(length(ni_surv)), function(i) matrix(Wtime[i, ], nrow = 1))

  # ----------------------------------------------------------
  # Initial values
  # ----------------------------------------------------------
  # Helper: ensure D is symmetric and strictly positive definite before any
  # arma::inv_sympd(D) call inside ll_lla_tv.cpp.  For univariate fits with
  # near-flat biomarkers (e.g., albumin), the LLA intercept/slope columns are
  # nearly collinear and the empirical Bsigma (or D.new from M-step) can be
  # PSD-but-singular -> inv_sympd() then aborts the EM at the first iteration.
  pd_safe <- function(M, eps = 1e-8) {
    M <- 0.5 * (M + t(M))
    ev <- tryCatch(eigen(M, symmetric = TRUE, only.values = TRUE)$values,
                   error = function(e) NA_real_)
    if (anyNA(ev) || min(ev) < eps) {
      shift <- if (anyNA(ev)) eps else (eps - min(ev) + eps)
      M <- M + shift * diag(nrow(M))
    }
    M
  }

  D     <- pd_safe(as.matrix(Matrix::bdiag(Bsigma)))
  var.e <- Ysigma2

  eta <- surv.init$eta
  phi <- surv.init$phi

  # ----------------------------------------------------------
  # Survival-related precomputation
  # ----------------------------------------------------------
  sv <- surv.mod(surv.init$ph, data, l0.init = surv.init$lamb)

  ft  <- sv$ft
  nev <- sv$nev

  Di <- sv$Di
  Deltai.list <- as.list(Di)

  l0  <- sv$l0
  l0i <- sv$l0i
  l0i.list <- as.list(l0i)
  l0u <- sv$l0u

  Fi <- sv$Fi
  surv.times <- sv$surv.times

  # ----------------------------------------------------------
  # B-spline basis pre-computation (from RefinedfastEM_tv)
  # ----------------------------------------------------------
  BB <- list()
  for (k in seq_len(n_LLA)) {
    Null_check <- unlist(B_list[[k]])
    if (is.null(Null_check)) {
      BB[[k]] <- matrix(rep(1, length(ft)), ncol = 1)
    } else {
      BB[[k]] <- predict(B_list[[k]], ft)
    }
  }

  eta.inds <- sapply(BB, function(x) 2 * ncol(x))

  BBc <- do.call(cbind, BB)
  BBcd <- do.call(cbind, lapply(seq_len(ncol(BBc)), function(i) {
    cbind(BBc[, i], BBc[, i])
  }))

  BBi <- lapply(surv.times, function(x) matrix(BBcd[x, ], ncol = length(eta_n)))

  # ----------------------------------------------------------
  # Working V
  # ----------------------------------------------------------
  V <- lapply(seq_len(n), function(ii) {
    mi_i <- mi[[ii]]
    blocks <- lapply(seq_len(nK), function(k) diag(var.e[k], nrow = mi_i[k], ncol = mi_i[k]))
    as.matrix(Matrix::bdiag(blocks))
  })

  # ----------------------------------------------------------
  # Convert LLA initial random effects to subject-wise vectors
  # (identical to RefinedfastEM_LLA)
  # ----------------------------------------------------------
  make_b_subject_list <- function(bLLA, nsub) {
    stopifnot(is.list(bLLA), length(bLLA) >= 1)
    KK <- length(bLLA)
    for (k in seq_len(KK)) {
      if (!is.matrix(bLLA[[k]]) || nrow(bLLA[[k]]) != nsub || ncol(bLLA[[k]]) != 2) {
        stop(sprintf("bLLA[[%d]] must be a %d x 2 matrix.", k, nsub))
      }
    }
    nm <- unlist(lapply(seq_len(KK), function(k) paste0("b", k, c("_0", "_1"))))
    b <- lapply(seq_len(nsub), function(i) {
      v <- as.numeric(unlist(lapply(seq_len(KK), function(k) bLLA[[k]][i, ])))
      names(v) <- nm
      v
    })
    b
  }

  make_c_subject_list <- function(cLLA, nsub) {
    stopifnot(is.list(cLLA), length(cLLA) >= 1)
    KK <- length(cLLA)
    c0 <- as.numeric(unlist(lapply(seq_len(KK), function(k) {
      ck <- cLLA[[k]]
      if (is.matrix(ck) || is.data.frame(ck)) {
        if (nrow(ck) != 1 || ncol(ck) != 2) stop(sprintf("cLLA[[%d]] must be 1x2.", k))
        as.numeric(ck[1, ])
      } else {
        ck <- as.numeric(ck)
        if (length(ck) != 2) stop(sprintf("cLLA[[%d]] must be length 2.", k))
        ck
      }
    })))
    nm <- unlist(lapply(seq_len(KK), function(k) paste0("c", k, c("_0", "_1"))))
    names(c0) <- nm
    replicate(nsub, c0, simplify = FALSE)
  }

  b <- make_b_subject_list(bLLA, nsub = n_subject)
  c <- make_c_subject_list(cLLA, nsub = n_subject)

  # ----------------------------------------------------------
  # Parameter vector bookkeeping
  # ----------------------------------------------------------
  b.inds <- split(seq_len(2 * nK), rep(seq_len(nK), each = 2))

  vD <- vech(D)
  names(vD) <- paste0("D[", apply(which(lower.tri(D, TRUE), arr.ind = TRUE), 1, paste0, collapse = ","), "]")

  params <- c(vD, var.e, eta, phi)

  eta_names <- eta_n  # already constructed by InitVal_LLAJEL_tv

  name_vec <- c(
    names(vD),
    paste0("var.e_", seq_len(nK)),
    eta_names
  )

  if (!is.null(phi) && length(phi) > 0) {
    name_vec <- c(name_vec, paste0("phi_", seq_along(phi)))
  }

  names(params) <- name_vec

  # Data bundle
  dmats <- list(
    Y = Y, Z = Z,
    W_kernel = W_kernel,
    mi = mi,
    Yk = Yk, Zk = Zk,
    w_longK = w_longK,
    K = K, Deltai = Deltai.list,
    BB = BB, BBi = BBi,
    c = c
  )

  if (collect.hist) {
    iter.hist <- data.frame(iter = iter, t(params))
  }

  # GH quadrature
  gh <- statmod::gauss.quad.prob(gh.nodes, "normal")
  v <- gh$nodes
  w <- gh$weights

  message("Starting EM Algorithm (LJM, time-varying)")

  EM.time <- numeric(0)

  # ----------------------------------------------------------
  # EM loop
  # ----------------------------------------------------------
  while (diff > tol && iter < max.iter) {

    p1 <- proc.time()[3]

    # ======================
    # E-step: posterior mode b.hat and covariance Sigmai
    # Uses ll_lla_tv / gradll_lla_tv / sdll_lla_tv
    # ======================

    # Compute etaBBs for current eta
    etaBBs <- etaBBs_compute(BBi, eta, eta.inds)

    b.hat <- mapply(
      function(b, Y, Z, V, K, Delta, l0i, l0u, wker, c, etaBB) {

        mi_i <- length(Y)

        ucminf::ucminf(
          par = b,
          fn  = ll_lla_tv,
          gr  = gradll_lla_tv,
          Y = Y, Z = Z, V = V, D = D,
          mi = mi_i,
          K = K,
          Delta = Delta,
          l0i = l0i,
          l0u = l0u,
          phi = phi,
          etaBB = etaBB,
          w = wker,
          c = c,
          control = list(xtol = 1e-3, grtol = 1e-6)
        )$par
      },
      b = b, Y = Y, Z = Z, V = V, K = K,
      Delta = Deltai.list,
      l0i = l0i.list,
      l0u = l0u,
      wker = W_kernel,
      c = c,
      etaBB = etaBBs,
      SIMPLIFY = FALSE
    )

    # ======================
    # c update (population center of the random-effects distribution)
    # ======================
    # Manuscript M-step: hat{c} = n_s^{-1} sum_i hat{b}_i.  update_c = FALSE
    # reproduces the legacy frozen-c behavior.  (Same as RefinedfastEM_LLA.)
    if (isTRUE(update_c)) {
      c_new_vec <- Reduce("+", b.hat) / length(b.hat)
      names(c_new_vec) <- names(c[[1]])
      c <- replicate(length(b.hat), c_new_vec, simplify = FALSE)
    }

    # Posterior covariance at mode
    Sigmai <- mapply(
      function(bhat, Z, V, K, l0u, wker, etaBB) {
        solve(-1 * sdll_lla_tv(bhat, Z, D, V, K, l0u, phi, etaBB, wker))
      },
      bhat = b.hat, Z = Z, V = V, K = K, l0u = l0u, wker = W_kernel, etaBB = etaBBs,
      SIMPLIFY = FALSE
    )

    S <- lapply(Sigmai, function(Si) lapply(b.inds, function(ix) Si[ix, ix]))

    D.newi <- mapply(function(Sigma_i, b_i, c_i) Sigma_i + tcrossprod(b_i - c_i),
                     Sigma_i = Sigmai, b_i = b.hat, c_i = c, SIMPLIFY = FALSE)

    # ======================
    # Kernel-weighted longitudinal expectations (unchanged from LLA)
    # ======================
    inv_omega_ker <- calc_inv_omega_kernel(mi, var.e, W_kernel)

    Sigma.longK <- calc_Sigma_longK(Zk, S, nK)
    mu.longK    <- calc_mu_longK_noX(Yk, Zk, b.hat, b.inds, nK)
    Ewe         <- calc_Ee_kernel(mu.longK, Sigma.longK, w_longK, nK)

    # ======================
    # Survival expectations (time-varying)
    # ======================
    mu_surv     <- calc_mu_surv_noX_tv(Y, Z, inv_omega_ker, K, D, phi, etaBBs, c)
    Sigma2_surv <- calc_Sigma2_surv_noX_tv(Z, inv_omega_ker, D, etaBBs)

    Es_exp <- Esurv_exp_t(w, v,
                          mu = mu_surv, variance = Sigma2_surv,
                          mu_new = mu_surv, variance_new = Sigma2_surv,
                          l0i = l0i, l0u = l0u)

    # Newton items for (eta, phi)
    Sge <- Setaphi_lla_tv(c(eta, phi),
                          Y, Z, inv_omega_ker,
                          K, D,
                          l0i, l0u, Di,
                          nK, w, v,
                          mu_surv, Sigma2_surv,
                          BBi, eta.inds,
                          eps = 0.0001,
                          c_list = c)

    Hge <- Hetaphi_lla_tv(c(eta, phi),
                          Y, Z, inv_omega_ker,
                          K, D,
                          l0i, l0u, Di,
                          nK, w, v,
                          mu_surv, Sigma2_surv,
                          BBi, eta.inds,
                          eps = 0.0001,
                          c_list = c)

    # ======================
    # M-step
    # ======================

    # D update
    D.new <- Reduce("+", D.newi) / n
    D.new <- pd_safe(D.new)   # keep symmetric + strictly PD for inv_sympd

    # var.e update (kernel-weighted)
    num <- colSums(do.call(rbind, Ewe))
    den <- colSums(do.call(rbind, lapply(w_longK, function(wi) sapply(wi, sum))))
    var.e.new <- num / den

    # Baseline hazard update (time-varying version)
    lambda <- update_lambda_tv(Es_exp, l0, n)
    l0.new <- nev / rowSums(do.call(cbind, lambda))
    l0u.new <- lapply(l0u, function(x) {
      ll <- length(x)
      l0.new[seq_len(ll)]
    })
    l0i.new <- numeric(length(Di))
    l0i.new[Di == 0] <- 0
    l0i.new[Di == 1] <- l0.new[match(Fi[Di == 1, 2], ft)]

    # (eta, phi) update (Newton step)
    eta.phi.new <- c(eta, phi) - solve(Hge, Sge)

    q_eta <- length(eta)

    eta.new <- eta.phi.new[seq_len(q_eta)]

    if (length(phi) > 0) {
      phi.new <- eta.phi.new[(q_eta + 1):length(eta.phi.new)]
    } else {
      phi.new <- numeric(0)
    }

    EM.time[iter + 1] <- proc.time()[3] - p1

    # ======================
    # Convergence diagnostics
    # ======================
    params.new <- c(vech(D.new), var.e.new, eta.new, phi.new)
    names(params.new) <- names(params)

    if (verbose) print(round(params.new, 4))

    if (diff.type == "abs") {
      diffs  <- abs(params.new - params)
      b.diff <- max(abs(do.call(rbind, b.hat) - do.call(rbind, b)))
    } else if (diff.type == "abs.rel") {
      diffs  <- abs(params.new - params) / (abs(params) + 1e-3)
      b.diff <- max(abs(do.call(rbind, b.hat) - do.call(rbind, b)) /
                      (abs(do.call(rbind, b)) + 1e-3))
    } else {
      stop("diff.type must be 'abs' or 'abs.rel'")
    }

    diff <- max(diffs)

    message("\nIteration ", iter + 1, " max diff: ", round(diff, 5))
    message("Largest change: ", names(params)[which.max(diffs)])
    message("--> old: ", params[which.max(diffs)], " new: ", params.new[which.max(diffs)])
    message("Largest change in random effects: ", round(b.diff, 3))

    # ======================
    # Update for next iteration
    # ======================
    params <- params.new
    D     <- D.new
    var.e <- var.e.new

    V <- lapply(mi, function(mi_i) {
      blocks <- lapply(seq_len(nK), function(k) diag(var.e[k], nrow = mi_i[k], ncol = mi_i[k]))
      as.matrix(Matrix::bdiag(blocks))
    })

    eta <- eta.new
    phi <- phi.new

    b <- b.hat

    l0  <- l0.new
    l0u <- l0u.new
    l0i <- l0i.new
    l0i.list <- as.list(l0i)

    iter <- iter + 1

    if (collect.hist) iter.hist <- rbind(iter.hist, c(iter = iter, t(params)))
  }

  # ----------------------------------------------------------
  # Final formatting
  # ----------------------------------------------------------
  rownames(D) <- sapply(seq_along(b.inds), function(k) paste0("b", k, "_", seq_along(b.inds[[k]]) - 1))
  colnames(D) <- rownames(D)

  names(var.e) <- paste0("var.e_", seq_len(nK))
  names(eta)   <- eta_n
  if (!is.null(phi) && length(phi) > 0L) {
    names(phi) <- paste0("phi_", seq_along(phi))
  } else {
    phi <- numeric(0)
  }

  coeffs <- list(
    beta   = NULL,
    var.e  = var.e,
    D      = D,
    eta    = eta,
    phi    = phi,
    c      = c[[1]],   # fitted center of the random-effects distribution
    hazard = cbind(ft, l0),
    l0u    = l0u
  )

  rtn <- list(
    REs   = do.call(rbind, b),
    coeffs = coeffs,
    dmats  = dmats,
    EMtime = round(sum(EM.time), 2),
    comp.time = round(proc.time()[3] - start.time, 2)
  )

  if (collect.hist) rtn$history <- iter.hist

  # ----------------------------------------------------------
  # Optional post-fit calculations: Vcov via PRES_hessian_LLA_tv
  # ----------------------------------------------------------
  if (post.process) {
    message("\nStarting post-fit calculations...")
    pp.start <- proc.time()[3]

    # Recompute etaBBs at final eta
    etaBBs <- etaBBs_compute(BBi, eta, eta.inds)

    b.hat <- mapply(
      function(b, Y, Z, V, K, Delta, l0i, l0u, wker, c, etaBB) {
        mi_i <- length(Y)
        ucminf::ucminf(
          par = b,
          fn  = ll_lla_tv,
          gr  = gradll_lla_tv,
          Y = Y, Z = Z, V = V, D = D,
          mi = mi_i,
          K = K,
          Delta = Delta,
          l0i = l0i,
          l0u = l0u,
          phi = phi,
          etaBB = etaBB,
          w = wker,
          c = c,
          control = list(xtol = 1e-3, grtol = 1e-6)
        )$par
      },
      b = b,
      Y = Y, Z = Z, V = V, K = K,
      Delta = Deltai.list, l0i = l0i.list, l0u = l0u,
      wker = W_kernel, c = c, etaBB = etaBBs,
      SIMPLIFY = FALSE
    )

    Sigmai <- mapply(
      function(bhat, Z, V, K, l0u, wker, etaBB) {
        solve(-1 * sdll_lla_tv(bhat, Z, D, V, K, l0u, phi, etaBB, wker))
      },
      bhat = b.hat, Z = Z, V = V, K = K, l0u = l0u, wker = W_kernel, etaBB = etaBBs,
      SIMPLIFY = FALSE
    )

    S <- lapply(Sigmai, function(Si) lapply(b.inds, function(ix) Si[ix, ix]))

    dmats$c <- c   # refresh with the fitted center (see RefinedfastEM_LLA)

    H <- PRES_hessian_LLA_tv(coeffs, dmats, V, b, b.hat, Sigmai, S,
                              l0i, l0u, gh.nodes, n, q, nK, eta.inds, nev, Fi, delta = 0.0001)

    # Regularized inverse with robust fallback
    if (any(!is.finite(H))) {
      warning("PRES Hessian contains NaN/Inf entries; replacing with 0.")
      H[!is.finite(H)] <- 0
    }

    # Eigen-based regularized inverse of the information matrix I = -H.
    # Floor the eigenvalues of I at an ABSOLUTE lambda (NOT tol*max|eig|).
    # Rationale: small variance components in D have legitimately HUGE Fisher
    # information (~1/sigma^4 -> eigenvalues ~1e12); a max-relative floor
    # (lam = tol*max|eig| ~ 1e7) then crushes the variance of EVERY direction
    # -- including the well-identified survival eta's (O(1) eigenvalues) --
    # to ~0, producing spuriously tiny SE.  An absolute floor instead:
    #   - well-identified directions (eta, large-info D)  -> 1/eigenvalue  (correct)
    #   - weakly/unidentified directions (eig < lam)      -> 1/lam = large (honest wide bands)
    # Eigen path is numerically stable even when cond(H) is huge (unlike solve()).
    # Tikhonov (diagonal) ridge on the information matrix: V = (I + lam*Identity)^{-1},
    # implemented via eigen for numerical stability (the legitimate huge D-variance
    # eigenvalue keeps cond(H) ~1e15, so a direct solve(-H + lam*I) would lose
    # precision). Indefinite directions are clipped to 0 before adding lam.
    #   well-identified dir (e_i >> lam) -> 1/e_i   (correct)
    #   null dir         (e_i ~ 0)       -> 1/lam   (finite, wide)
    ridge_vcov <- function(H, lam = 1e-6) {
      I  <- -H; I <- 0.5 * (I + t(I))
      ev <- eigen(I, symmetric = TRUE)
      Vr <- ev$vectors %*% diag(1 / (pmax(ev$values, 0) + lam)) %*% t(ev$vectors)
      0.5 * (Vr + t(Vr))
    }
    Vcov <- ridge_vcov(H)

    rownames(Vcov) <- names(params)
    colnames(Vcov) <- names(params)
    rtn$Hessian <- H

    pp.end <- proc.time()[3]
    message("\nDone")

    rtn$Vcov <- Vcov
    rtn$REs  <- do.call(rbind, b)
    rtn$postprocess.time <- round(pp.end - pp.start, 2)
    rtn$comp.time <- round(proc.time()[3] - start.time, 2)
  }

  rtn
}
