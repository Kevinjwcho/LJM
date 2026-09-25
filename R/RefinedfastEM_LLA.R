# ============================================================
# Refined Fast EM for Kernel-weighted LJM (Fixed effects)
# ------------------------------------------------------------
# - Longitudinal part: Kernel-weighted Local Linear Approximation (LLA)
#   * No X / No beta in the longitudinal likelihood
#   * Kernel weights enter only through the weighted residual likelihood
# - Survival part: NOT kernel-weighted
#   * Kernel affects survival only through the LOCAL posterior (b_i | Y_i)
# - Updates: D, var.e, baseline hazard l0, (eta, phi)
# - Optional post-processing: observed information / Vcov via PRES_hessian
#
# Required external functions (already in your project):
#   surv.mod(), update_lambda(), vech()
#   ll_lla(), gradll_lla(), sdll_lla()
#   calc_inv_omega_kernel(), calc_Sigma_longK(), calc_mu_longK_noX(), calc_Ee_kernel()
#   calc_mu_surv_noX(), calc_Sigma2_surv_noX()
#   Esurv_exp(), Setaphi(), Hetaphi()
#   score_JEL(), PRES_hessian()  (your rewritten versions: no X/beta, kernel posterior)
# ============================================================

RefinedfastEM_LLA <- function(
    data,
    n_LLA, ID_surv, ni_surv,
    Z.st, Y.st, Wker.st, Wtime,
    Ysigma2, Bsigma,
    bLLA, cLLA, n_subject,
    surv.init, ni, ncw,
    gh.nodes = 3, collect.hist = TRUE, max.iter = 200,
    tol = 0.01, diff.type = "abs.rel",
    post.process = FALSE, verbose = FALSE,
    update_c = TRUE,
    rho = 0          # transformation parameter of G (0 = Cox)
){
  
  # ----------------------------------------------------------
  # Basic setup
  # ----------------------------------------------------------
  start.time <- proc.time()[3]
  
  nK <- n_LLA                 # number of longitudinal processes
  q  <- 2 * nK                # random effects dimension (intercept+slope per process)
  uids <- ID_surv
  n <- length(uids)
  
  diff   <- 100
  b.diff <- 100
  iter   <- 0
  
  # ----------------------------------------------------------
  # Construct subject-wise stacked longitudinal objects
  #   Z[[i]]: block-diagonal design matrix (stacked over k)
  #   Y[[i]]: stacked response vector (stacked over k)
  #   mi[[i]]: vector of observation counts (m_{ik})
  #   Zk[[i]][[k]], Yk[[i]][,k]: per-process objects
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
  
  # Kernel weights:
  # - W_kernel[[i]]: stacked weights aligned with stacked Y[[i]]
  # - w_longK[[i]][[k]]: weights for process k aligned with Yk[[i]][,k]
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
  
  # Survival design rows (one row per subject; used in K*phi)
  K <- lapply(seq_len(length(ni_surv)), function(i) matrix(Wtime[i, ], nrow = 1))
  
  # ----------------------------------------------------------
  # Initial values
  # ----------------------------------------------------------
  # Helper: ensure D is symmetric and strictly positive definite before any
  # arma::inv_sympd(D) call inside ll_lla.cpp.  See RefinedfastEM_LLA_tv.R for
  # the rationale (univariate near-collinear b columns -> singular D).
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

  D     <- pd_safe(as.matrix(Matrix::bdiag(Bsigma)))  # random effects covariance
  var.e <- Ysigma2                           # residual variances (length nK)
  
  # Survival parameters (phi, eta)
  eta <- surv.init$eta
  phi <- surv.init$phi
  
  # ----------------------------------------------------------
  # Survival-related precomputation
  # ----------------------------------------------------------
  sv <- surv.mod(surv.init$ph, data, l0.init = surv.init$lamb)
  
  ft <- sv$ft
  nev <- sv$nev
  
  Di <- sv$Di
  Deltai.list <- as.list(Di)
  
  l0  <- sv$l0
  l0i <- sv$l0i
  l0i.list <- as.list(l0i)
  l0u <- sv$l0u
  
  Fi <- sv$Fi
  
  # ----------------------------------------------------------
  # Working V (block-diagonal var.e replicated by mi)
  # Note: ll_lla may not need V if it uses inv_omega_ker internally,
  # but you already pass V consistently; keep it.
  #
  # Each block is diag(var.e[k], mi_i[k]), so the resulting
  # block-diagonal is itself diagonal.  diag(rep(...)) gives a
  # bit-identical result with no sparse-coercion overhead
  # (replaces Matrix::bdiag() which triggers isSymmetric /
  #  all.equal chains that dominated the K=5 profile).
  # ----------------------------------------------------------
  V <- lapply(mi, function(mi_i) diag(rep(var.e, times = mi_i)))
  
  # ----------------------------------------------------------
  # Convert LLA initial random effects to subject-wise vectors
  # bLLA: list length nK, each is n x 2 matrix
  # cLLA: list length nK, each is 1x2 (global center)
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
  
  # NOTE: these objects must exist in your environment:
  #   bLLA (list length nK) and cLLA (list length nK)
  #   n_subject
  b <- make_b_subject_list(bLLA, nsub = n_subject)
  c <- make_c_subject_list(cLLA, nsub = n_subject)
  
  # ----------------------------------------------------------
  # Parameter vector bookkeeping
  # ----------------------------------------------------------
  b.inds <- split(seq_len(2 * nK), rep(seq_len(nK), each = 2))
  
  vD <- vech(D)
  names(vD) <- paste0("D[", apply(which(lower.tri(D, TRUE), arr.ind = TRUE), 1, paste0, collapse = ","), "]")
  
  params <- c(vD, var.e, eta, phi)
  # names(params) <- c(
  #   names(vD),
  #   paste0("var.e_", seq_len(nK)),
  #   sapply(seq_along(b.inds), function(k) paste0("eta_", k, "_", seq_along(b.inds[[k]]))),
  #   paste0("phi_", seq_along(phi))
  # )
  
  eta_names <- unlist(
    lapply(seq_along(b.inds), function(k) {
      paste0("eta_", k, "_", seq_along(b.inds[[k]]))
    }),
    use.names = FALSE
  )
  
  name_vec <- c(
    names(vD),
    paste0("var.e_", seq_len(nK)),
    eta_names
  )
  
  if (!is.null(phi) && length(phi) > 0) {
    name_vec <- c(name_vec, paste0("phi_", seq_along(phi)))
  }
  
  names(params) <- name_vec
  
  # Data bundle used in PRES score/Hessian
  dmats <- list(
    Y = Y, Z = Z,
    W_kernel = W_kernel,
    mi = mi,
    Yk = Yk, Zk = Zk,
    w_longK = w_longK,
    K = K, Deltai = Deltai.list,
    c = c
  )
  
  if (collect.hist) {
    iter.hist <- data.frame(iter = iter, t(params))
  }
  
  # GH quadrature nodes/weights
  gh <- statmod::gauss.quad.prob(gh.nodes, "normal")
  v <- gh$nodes
  w <- gh$weights
  
  message("Starting EM Algorithm (LJM)")

  EM.time <- numeric(0)

  # per-block wall-time accumulators (diagnostic: locate K-scaling bottleneck)
  prof <- list(
    bhat         = 0,   # E-step: ucminf mode-finding per subject
    Sigmai       = 0,   # E-step: solve(-sdll_lla) per subject + S split + D.newi
    expectations = 0,   # E-step: calc_* (long + surv) + Esurv_exp
    Setaphi      = 0,   # E-step: Rcpp score for (eta,phi)
    Hetaphi      = 0,   # E-step: Rcpp Hessian for (eta,phi)
    Mstep        = 0    # M-step: D, var.e, lambda, eta/phi Newton
  )

  # ----------------------------------------------------------
  # EM loop
  # ----------------------------------------------------------
  guard_hits <- 0L; guard_reset_done <- FALSE   # (eta, phi) Newton guard, see below
  while (diff > tol && iter < max.iter) {
    
    p1 <- proc.time()[3]
    
    # ======================
    # E-step: posterior mode b.hat and covariance Sigmai
    # ======================
    .t <- proc.time()[3]
    b.hat <- mapply(
      function(b, Y, Z, V, K, Delta, l0i, l0u, wker, c) {

        mi_i <- length(Y)  # stacked length across all k

        ucminf::ucminf(
          par = b,
          fn  = ll_lla,
          gr  = gradll_lla,
          Y = Y, Z = Z, V = V, D = D,
          mi = mi_i,
          K = K,
          Delta = Delta,
          l0i = l0i,
          l0u = l0u,
          phi = phi,
          eta = eta,
          w = wker,
          c = c,
          rho = rho,
          control = list(xtol = 1e-3, grtol = 1e-6)
        )$par
      },
      b = b, Y = Y, Z = Z, V = V, K = K,
      Delta = Deltai.list,
      l0i = l0i.list,
      l0u = l0u,
      wker = W_kernel,
      c = c,
      SIMPLIFY = FALSE
    )
    prof$bhat <- prof$bhat + (proc.time()[3] - .t)

    # ======================
    # c update (population center of the random-effects distribution)
    # ======================
    # Manuscript M-step: hat{c} = n_s^{-1} sum_i hat{b}_i, evaluated at the
    # posterior means hat{b}_i. This was previously frozen at the kernel-LLA
    # initial center (cLLA), which biased both c and the D update (D is formed
    # from deviations b_i - c). update_c = FALSE reproduces the legacy behavior.
    if (isTRUE(update_c)) {
      c_new_vec <- Reduce("+", b.hat) / length(b.hat)
      names(c_new_vec) <- names(c[[1]])
      c <- replicate(length(b.hat), c_new_vec, simplify = FALSE)
    }

    # Posterior covariance approximation at mode
    .t <- proc.time()[3]
    Sigmai <- mapply(
      function(bhat, Z, V, K, l0u, wker, Delta) {
        solve(-1 * sdll_lla(bhat, Z, D, V, K, l0u, phi, eta, wker, Delta = Delta, rho = rho))
      },
      bhat = b.hat, Z = Z, V = V, K = K, l0u = l0u, wker = W_kernel, Delta = Deltai.list,
      SIMPLIFY = FALSE
    )

    # Split Sigmai into process-specific 2x2 blocks (needed for longitudinal Ewe)
    S <- lapply(Sigmai, function(Si) lapply(b.inds, function(ix) Si[ix, ix]))

    # E[b b'] piece for updating D
    D.newi <- mapply(function(Sigma_i, b_i, c_i) Sigma_i + tcrossprod(b_i - c_i),
                     Sigma_i = Sigmai, b_i = b.hat, c_i = c, SIMPLIFY = FALSE)
    prof$Sigmai <- prof$Sigmai + (proc.time()[3] - .t)
    
    # ======================
    # Kernel-weighted longitudinal expectations for var.e
    # ======================
    .t <- proc.time()[3]
    inv_omega_ker <- calc_inv_omega_kernel(mi, var.e, W_kernel)

    Sigma.longK <- calc_Sigma_longK(Zk, S, nK)
    mu.longK    <- calc_mu_longK_noX(Yk, Zk, b.hat, b.inds, nK)
    Ewe         <- calc_Ee_kernel(mu.longK, Sigma.longK, w_longK, nK)

    # ======================
    # Survival expectations using LOCAL posterior (kernel only via inv_omega_ker)
    # ======================
    mu_surv     <- calc_mu_surv_noX(Y, Z, inv_omega_ker, K, D, phi, eta, c)
    Sigma2_surv <- calc_Sigma2_surv_noX(Z, inv_omega_ker, D, eta)

    Es_exp <- Esurv_exp(w, v,
                        mu = mu_surv, variance = Sigma2_surv,
                        mu_new = mu_surv, variance_new = Sigma2_surv,
                        l0i = l0i, l0u = l0u, rho = rho)
    prof$expectations <- prof$expectations + (proc.time()[3] - .t)

    # Newton items for (eta, phi)
    .t <- proc.time()[3]
    Sge <- Setaphi(c(eta, phi),
                   Y, Z, inv_omega_ker,
                   K, D,
                   l0i, l0u, Di,
                   nK, w, v,
                   mu_surv, Sigma2_surv,
                   eps = 0.0001,
                   c_list = c, rho = rho)
    prof$Setaphi <- prof$Setaphi + (proc.time()[3] - .t)

    .t <- proc.time()[3]
    Hge <- Hetaphi(c(eta, phi),
                   Y, Z, inv_omega_ker,
                   K, D,
                   l0i, l0u, Di,
                   nK, w, v,
                   mu_surv, Sigma2_surv,
                   eps = 0.0001,
                   c_list = c, rho = rho)
    prof$Hetaphi <- prof$Hetaphi + (proc.time()[3] - .t)

    if (verbose && (!all(is.finite(Sge)) || !all(is.finite(Hge)))) {
      .fin <- function(x) all(is.finite(unlist(x)))
      .rng <- function(x) { x <- unlist(x); x <- x[is.finite(x)]; if (length(x)) sprintf("[%.3g, %.3g]", min(x), max(x)) else "none finite" }
      message(sprintf("  [guard diag] finite: b.hat=%s Sigmai=%s D.newi=%s Ewe=%s mu_surv=%s Sigma2_surv=%s Es_exp=%s l0u=%s l0i=%s | ranges: b.hat %s, mu_surv %s, Sigma2_surv %s, Es_exp %s, l0u %s",
        .fin(b.hat), .fin(Sigmai), .fin(D.newi), .fin(Ewe), .fin(mu_surv), .fin(Sigma2_surv), .fin(Es_exp), .fin(l0u), .fin(l0i),
        .rng(b.hat), .rng(mu_surv), .rng(Sigma2_surv), .rng(Es_exp), .rng(l0u)))
    }
    
    # ======================
    # M-step
    # ======================
    .t <- proc.time()[3]

    # D update
    D.new <- Reduce("+", D.newi) / n
    D.new <- pd_safe(D.new)   # keep symmetric + strictly PD for inv_sympd

    # var.e update: numerator = sum_i Ewe_{ik}, denominator = sum_i sum_r w_{ikr}
    num <- colSums(do.call(rbind, Ewe))
    den <- colSums(do.call(rbind, lapply(w_longK, function(wi) sapply(wi, sum))))
    var.e.new <- num / den

    # Baseline hazard update
    lambda <- update_lambda(Es_exp, l0u, l0, n)
    l0.new <- nev / rowSums(do.call(cbind, lambda))
    l0u.new <- lapply(l0u, function(x) {
      ll <- length(x)
      l0.new[seq_len(ll)]
    })
    l0i.new <- numeric(length(Di))
    l0i.new[Di == 0] <- 0
    l0i.new[Di == 1] <- l0.new[match(Fi[Di == 1, 2], ft)]

    # (eta, phi) update (Newton step)
    # eta.phi.new <- c(eta, phi) - solve(Hge, Sge)
    # eta.new <- eta.phi.new[1:(2 * nK)]
    # phi.new <- eta.phi.new[(2 * nK + 1):length(eta.phi.new)]

    # Guarded Newton step (JEL 2.3). On small risk sets or badly scaled markers
    # the Cox-based initial (eta, phi) can be so extreme that the score and
    # Hessian are non-finite (exp overflow), which turned (eta, phi) into NaN and
    # stopped the EM with an uninformative error.
    # (i) non-finite score/Hessian or singular Hessian: the first time this
    #     happens (eta, phi) is restarted from zero, where the score is finite;
    #     afterwards the update is skipped for that iteration;
    # (ii) otherwise the largest coordinate move is capped at `max_step`.
    # `guard_hits` counts interventions and is returned in the fit.
    max_step <- 10
    nr_step <- if (all(is.finite(Sge)) && all(is.finite(Hge)))
      tryCatch(solve(Hge, Sge), error = function(e) rep(NA_real_, length(Sge)))
    else rep(NA_real_, length(Sge))
    if (!all(is.finite(nr_step))) {
      guard_hits <- guard_hits + 1L
      if (!guard_reset_done) {
        guard_reset_done <- TRUE
        if (verbose) message(sprintf(
          "Iteration %d: (eta, phi) Newton step not finite; restarting (eta, phi) from zero. Initial eta = %s",
          iter + 1, paste(signif(eta, 3), collapse = " ")))
        nr_step <- c(eta, phi)                     # theta - theta = 0
      } else {
        if (verbose) message(sprintf("Iteration %d: (eta, phi) Newton step not finite; update skipped", iter + 1))
        nr_step <- rep(0, length(Sge))
      }
    } else if (max(abs(nr_step)) > max_step) {
      guard_hits <- guard_hits + 1L
      if (verbose) message(sprintf("Iteration %d: (eta, phi) Newton step damped by %.3g",
                                   iter + 1, max_step / max(abs(nr_step))))
      nr_step <- nr_step * (max_step / max(abs(nr_step)))
    }
    eta.phi.new <- c(eta, phi) - nr_step

    q_eta <- length(eta)

    eta.new <- eta.phi.new[seq_len(q_eta)]

    if (length(phi) > 0) {
      phi.new <- eta.phi.new[(q_eta + 1):length(eta.phi.new)]
    } else {
      phi.new <- numeric(0)
    }
    prof$Mstep <- prof$Mstep + (proc.time()[3] - .t)

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
    
    # Update V (block-diagonal with var.e repeated by mi).
    # Diagonal shortcut: bit-identical to the old bdiag-of-diags
    # construction but skips the Matrix::bdiag -> sparse coerce ->
    # isSymmetric -> all.equal chain that was dominating per-iter cost.
    V <- lapply(mi, function(mi_i) diag(rep(var.e, times = mi_i)))
    
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
  names(eta)   <- sapply(seq_along(b.inds), function(k) paste0("eta", k, "_", seq_along(b.inds[[k]]) - 1))
  # names(phi)   <- paste0("phi_", seq_along(phi))
  if (!is.null(phi) && length(phi) > 0L) {
    names(phi) <- paste0("phi_", seq_along(phi))
  } else {
    phi <- numeric(0)
  }
  
  coeffs <- list(
    # beta kept for compatibility with downstream code, but not used in this model
    beta   = beta,
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
    EMtime = round(sum(EM.time), 2),
    prof = lapply(prof, function(x) round(x, 3)),
    EM.time = EM.time,
    comp.time = round(proc.time()[3] - start.time, 2)
  )
  
  if (collect.hist) rtn$history <- iter.hist
  
  # ----------------------------------------------------------
  # Optional post-fit calculations: Vcov via PRES_hessian
  # ----------------------------------------------------------
  if (post.process) {
    message("\nStarting post-fit calculations...")
    pp.start <- proc.time()[3]

    # dmats was built with the initial center; refresh with the fitted c so the
    # PRES Hessian / Vcov below are evaluated at the converged parameters.
    dmats$c <- c
    
    # Recompute b.hat at final parameters (same as E-step)
    b.hat <- mapply(
      function(b, Y, Z, V, K, Delta, l0i, l0u, wker, c) {
        
        mi_i <- length(Y)
        
        ucminf::ucminf(
          par = b,
          fn  = ll_lla,
          gr  = gradll_lla,
          Y = Y, Z = Z, V = V, D = D,
          mi = mi_i,
          K = K,
          Delta = Delta,
          l0i = l0i,
          l0u = l0u,
          phi = phi,
          eta = eta,
          w = wker,
          c = c,
          rho = rho,
          control = list(xtol = 1e-3, grtol = 1e-6)
        )$par
      },
      b = b,
      Y = Y, Z = Z, V = V, K = K,
      Delta = Deltai.list, l0i = l0i.list, l0u = l0u,
      wker = W_kernel, c = c,
      SIMPLIFY = FALSE
    )
    
    Sigmai <- mapply(
      function(bhat, Z, V, K, l0u, wker, Delta) {
        solve(-1 * sdll_lla(bhat, Z, D, V, K, l0u, phi, eta, wker, Delta = Delta, rho = rho))
      },
      bhat = b.hat, Z = Z, V = V, K = K, l0u = l0u, wker = W_kernel, Delta = Deltai.list,
      SIMPLIFY = FALSE
    )
    
    # IMPORTANT: dmats already includes W_kernel and w_longK
    H <- PRES_hessian(coeffs, dmats, V, b, b.hat, Sigmai, S,
                      l0i, l0u, gh.nodes, n, q, nK, nev, Fi, delta = 0.0001, rho = rho)

    # Handle NaN/Inf in Hessian
    if (any(!is.finite(H))) {
      warning("PRES Hessian contains NaN/Inf entries; replacing with 0.")
      H[!is.finite(H)] <- 0
    }

    # Eigen-based Tikhonov-ridged inverse of the information I = -H (ported
    # from the tv path, RefinedfastEM_LLA_tv.R). A plain solve(H) is singular
    # for smooth kernels; the old fallback -ginv(H) then returns spuriously ~0
    # SE. Flooring the information's eigenvalues at an ABSOLUTE lam instead:
    #   well-identified dir (e_i >> lam) -> 1/e_i   (correct)
    #   null/weak dir       (e_i ~ 0)    -> 1/lam   (finite, honest wide bands)
    ridge_vcov <- function(H, lam = 1e-6) {
      I  <- -H; I <- 0.5 * (I + t(I))
      ev <- eigen(I, symmetric = TRUE)
      Vr <- ev$vectors %*% diag(1 / (pmax(ev$values, 0) + lam)) %*% t(ev$vectors)
      0.5 * (Vr + t(Vr))
    }
    Vcov <- ridge_vcov(H)

    rownames(Vcov) <- names(params)
    colnames(Vcov) <- names(params)
    
    pp.end <- proc.time()[3]
    message("\nDone")
    
    rtn$Vcov <- Vcov
    Hn <- H; rownames(Hn) <- colnames(Hn) <- names(params)
    rtn$Hessian <- Hn      # raw PRES Hessian (diagnostic: compare Vcov methods)
    rtn$REs  <- do.call(rbind, b)
    rtn$postprocess.time <- round(pp.end - pp.start, 2)
    rtn$comp.time <- round(proc.time()[3] - start.time, 2)
  }

  rtn$newton_guard <- guard_hits   # number of guarded/damped (eta, phi) Newton steps (0 = untouched)
  rtn
}