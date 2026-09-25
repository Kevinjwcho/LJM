InitVal_LLAJEL <- function(
    # ---- OPTION A: bBLUP already given (list) ----
    bBLUP = NULL,              # list of matrices, each: n_subject x p_k
    
    # ---- OPTION B: if bBLUP is NULL, compute it via kernel LLA ----
    Y = NULL, Z = NULL, ID = NULL, W_kernel = NULL,  # lists
    n_subject = NULL, n_LME = NULL,
    ridge = 1e-8,
    
    # ---- survival side ----
    Indcs,
    start, stop, event,
    W, ncw,
    Wtime2,
    ID_surv
) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required.")
  }
  
  # -------------------------
  # 0) build/validate bBLUP
  # -------------------------
  if (is.null(bBLUP)) {
    if (is.null(Y) || is.null(Z) || is.null(ID) || is.null(W_kernel)) {
      stop("If bBLUP is NULL, you must provide Y, Z, ID, W_kernel.")
    }
    if (is.null(n_subject) || is.null(n_LME)) {
      stop("If bBLUP is NULL, you must provide n_subject and n_LME.")
    }
    if (length(Y) != n_LME || length(Z) != n_LME || length(ID) != n_LME || length(W_kernel) != n_LME) {
      stop("Y/Z/ID/W_kernel must be lists of length n_LME.")
    }
    
    # kernel local linear WLS per subject -> bBLUP[[k]]
    bBLUP <- vector("list", n_LME)
    for (k in seq_len(n_LME)) {
      yk  <- Y[[k]]
      Zk  <- Z[[k]]
      idk <- ID[[k]]
      wk  <- W_kernel[[k]]
      
      if (!is.matrix(Zk)) stop(sprintf("Z[[%d]] must be a matrix.", k))
      if (length(yk) != nrow(Zk) || length(idk) != nrow(Zk) || length(wk) != nrow(Zk)) {
        stop(sprintf("Length mismatch in (Y,Z,ID,W_kernel) for k=%d.", k))
      }
      
      pk <- ncol(Zk)
      bk <- matrix(0, nrow = n_subject, ncol = pk)
      
      for (i in seq_len(n_subject)) {
        idx <- which(idk == i & is.finite(wk) & wk > 0)
        if (length(idx) < pk) next
        
        Zi <- Zk[idx, , drop = FALSE]
        yi <- yk[idx]
        wi <- wk[idx]
        
        XtX <- crossprod(Zi, wi * Zi) + ridge * diag(pk)
        Xty <- crossprod(Zi, wi * yi)
        bk[i, ] <- as.numeric(solve(XtX, Xty))
      }
      bBLUP[[k]] <- bk
    }
  } else {
    if (!is.list(bBLUP) || length(bBLUP) < 1) stop("bBLUP must be a non-empty list.")
    if (is.null(n_LME)) n_LME <- length(bBLUP)
    if (length(bBLUP) != n_LME) stop("length(bBLUP) must equal n_LME.")
    if (is.null(n_subject)) n_subject <- nrow(bBLUP[[1]])
  }
  
  # -------------------------
  # 1) risk set unpack
  # -------------------------
  Index  <- Indcs$Index
  Index1 <- Indcs$Index1
  Index2 <- Indcs$Index2
  
  # -------------------------
  # 2) build Cox design (random only)
  #    time  : nSurv x n_eta
  #    time2 : M x n_eta
  # -------------------------
  n_eta <- sum(vapply(bBLUP, ncol, 1L))
  
  B_time  <- do.call(cbind, lapply(seq_len(n_LME), function(k) bBLUP[[k]][ID_surv, , drop = FALSE]))
  B_time2 <- do.call(cbind, lapply(seq_len(n_LME), function(k) bBLUP[[k]][Index,   , drop = FALSE]))
  
  eta_n <- c()
  for (k in seq_len(n_LME)) {
    pk <- ncol(bBLUP[[k]])
    eta_n <- c(eta_n, paste0("eta_b", k, "_", 0:(pk - 1)))
  }
  colnames(B_time)  <- eta_n
  colnames(B_time2) <- eta_n
  
  # -------------------------
  # 3) Cox fit: Surv ~ W + B_time
  # -------------------------
  data.init <- data.frame(start = start, stop = stop, event = event)
  
  if (ncw > 0) {
    if (!is.matrix(W) || ncol(W) != ncw) stop("W must be an nSurv x ncw matrix when ncw > 0.")
    Wdf <- as.data.frame(W)
    names(Wdf) <- paste0("W", seq_len(ncw))
    data.init <- cbind(data.init, Wdf)
  }
  
  data.init <- cbind(data.init, as.data.frame(B_time))
  
  rhs_terms <- c(if (ncw > 0) names(Wdf) else NULL, eta_n)
  fml <- stats::as.formula(
    paste0("survival::Surv(start, stop, event) ~ ", paste(rhs_terms, collapse = " + "))
  )
  
  fit <- survival::coxph(fml, data = data.init)
  
  phi.init <- if (ncw > 0) as.vector(fit$coefficients[seq_len(ncw)]) else numeric(0)
  eta.init <- as.vector(fit$coefficients[(ncw + 1):(ncw + n_eta)])

  # Safeguard (JEL 2.3): the unpenalised Cox fit on the initial local intercepts
  # and slopes can return extreme coefficients on small risk sets or badly scaled
  # markers (e.g. a slope coefficient of 20 on a raw-scale biomarker). The EM then
  # evaluates exp() of a linear predictor spanning hundreds of units, the first
  # E-step returns non-finite expectations and the fit fails. If the centred
  # initial linear predictor exceeds `lp_max` in absolute value, (eta, phi) are
  # shrunk proportionally; the EM moves them freely afterwards. The factor is
  # returned as `init_shrink` (1 = no intervention).
  lp_max <- 10
  init_shrink <- 1
  lp0 <- as.vector(B_time %*% eta.init)
  if (ncw > 0) lp0 <- lp0 + as.vector(W %*% phi.init)
  lp0 <- lp0[is.finite(lp0)]
  if (length(lp0)) {
    spread <- max(abs(lp0 - stats::median(lp0)))
    if (!is.finite(spread) || spread > lp_max) {
      init_shrink <- if (is.finite(spread)) lp_max / spread else 0
      message(sprintf("InitVal_LLAJEL: initial Cox (eta, phi) shrunk by %.3g (linear-predictor spread %.3g > %g)",
                      init_shrink, spread, lp_max))
      eta.init <- eta.init * init_shrink
      phi.init <- phi.init * init_shrink
    }
  }
  if (any(!is.finite(eta.init))) { eta.init[!is.finite(eta.init)] <- 0; init_shrink <- 0 }
  if (any(!is.finite(phi.init))) { phi.init[!is.finite(phi.init)] <- 0; init_shrink <- 0 }
  
  # -------------------------
  # 4) baseline hazard increments lambda
  # -------------------------
  M <- nrow(Wtime2)
  
  Wtime2_phi <- if (ncw > 0) as.vector(Wtime2 %*% phi.init) else rep(0, M)
  linpred <- if (n_eta == 1) as.vector(B_time2) * eta.init else as.vector(B_time2 %*% eta.init)
  
  temp <- exp(Wtime2_phi + linpred)
  
  lamb.init <- Index2 / calc_tapply_vect_sum(
    v1 = temp,
    v2 = as.integer(Index1 - 1)
  )
  
  list(
    phi   = phi.init,
    eta   = eta.init,
    lamb  = lamb.init,
    bBLUP = bBLUP,
    eta_n = eta_n,
    ph = fit,
    init_shrink = init_shrink
  )
}