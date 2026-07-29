
#=============== Initial Value Calculation for LLAJEL Time-Varying ===============#
# Extends InitVal_LLAJEL with B-spline time-varying association (eta)
# Mirrors InitVal_multiJEL_tv logic but uses bBLUP from LLA (not LME)

InitVal_LLAJEL_tv <- function(
    # ---- OPTION A: bBLUP already given (list) ----
    bBLUP = NULL,              # list of matrices, each: n_subject x p_k

    # ---- OPTION B: if bBLUP is NULL, compute it via kernel LLA ----
    Y = NULL, Z = NULL, ID = NULL, W_kernel = NULL,
    n_subject = NULL, n_LME = NULL,
    ridge = 1e-8,

    # ---- survival side ----
    Indcs,
    start, stop, event,
    W, ncw,
    Wtime2,
    ID_surv,

    # ---- B-spline specification ----
    Bs = NULL        # list of length n_LME, each: list(df, degree, knots, Bknots) or list(NULL)
) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required.")
  }

  # -------------------------
  # 0) build/validate bBLUP (same as InitVal_LLAJEL)
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
  # 2) B-spline construct (from InitVal_multiJEL_tv)
  # -------------------------
  if (is.null(Bs) || length(Bs) != n_LME) {
    stop("Bs must be a list of length n_LME for time-varying model.")
  }

  B_list <- list()
  for (j in seq_len(n_LME)) {
    Bs_config <- Bs[[j]]
    Null_check <- unlist(Bs_config)
    if (is.null(Null_check)) {
      B_list[[j]] <- list(NULL)
    } else {
      B_list[[j]] <- splines::bs(stop, df = Bs_config$df, knots = Bs_config$knots,
                                  degree = Bs_config$degree, intercept = TRUE) %>% as.matrix
    }
  }

  # Total number of eta parameters
  n_eta <- sapply(seq_len(n_LME), function(x) {
    Null_check <- unlist(B_list[[x]])
    if (is.null(Null_check)) {
      return(ncol(bBLUP[[x]]))
    } else {
      return(ncol(B_list[[x]]) * ncol(bBLUP[[x]]))
    }
  }) %>% sum

  # Build design matrices for Cox init
  fixedOrRand_mat  <- matrix(nrow = length(ID_surv), ncol = n_eta)
  fixedOrRand_mat2 <- matrix(nrow = length(Index), ncol = n_eta)

  # eta names (from InitVal_multiJEL_tv)
  eta_n <- sapply(seq_len(n_LME), function(x) {
    Null_check <- unlist(B_list[[x]])
    if (is.null(Null_check)) {
      paste0("eta_b", x, "_", 0:(ncol(bBLUP[[x]]) - 1))
    } else {
      outer(paste0("eta_B", seq_len(ncol(B_list[[x]]))),
            paste0("_b", x, "_", 0:(ncol(bBLUP[[x]]) - 1)),
            paste0) %>% t() %>% as.vector
    }
  }) %>% unlist

  var_loci <- 0
  for (k in seq_len(n_LME)) {
    Null_check <- unlist(B_list[[k]])
    if (is.null(Null_check)) {
      for (j in seq_len(ncol(bBLUP[[k]]))) {
        var_loci <- var_loci + 1
        fixedOrRand_mat[, var_loci]  <- bBLUP[[k]][ID_surv, j]
        fixedOrRand_mat2[, var_loci] <- bBLUP[[k]][Index, j]
      }
    } else {
      for (i in seq_len(ncol(B_list[[k]]))) {
        for (j in seq_len(ncol(bBLUP[[k]]))) {
          var_loci <- var_loci + 1
          fixedOrRand_mat[, var_loci]  <- B_list[[k]][, i] * bBLUP[[k]][ID_surv, j]
          fixedOrRand_mat2[, var_loci] <- B_list[[k]][Index, i] * bBLUP[[k]][Index, j]
        }
      }
    }
  }

  fixedOrRand.time  <- fixedOrRand_mat
  fixedOrRand.time2 <- fixedOrRand_mat2

  colnames(fixedOrRand.time)  <- eta_n
  colnames(fixedOrRand.time2) <- eta_n

  # -------------------------
  # 3) Cox fit: Surv ~ W + fixedOrRand.time
  # -------------------------
  data.init <- data.frame(start = start, stop = stop, event = event)

  if (ncw > 0) {
    if (!is.matrix(W) || ncol(W) != ncw) stop("W must be an nSurv x ncw matrix when ncw > 0.")
    Wdf <- as.data.frame(W)
    names(Wdf) <- paste0("W", seq_len(ncw))
    data.init <- cbind(data.init, Wdf)
  }

  data.init <- cbind(data.init, as.data.frame(fixedOrRand.time))

  rhs_terms <- c(if (ncw > 0) names(Wdf) else NULL, eta_n)
  fml <- stats::as.formula(
    paste0("survival::Surv(start, stop, event) ~ ", paste(rhs_terms, collapse = " + "))
  )

  fit <- survival::coxph(fml, data = data.init)

  phi.init <- if (ncw > 0) as.vector(fit$coefficients[seq_len(ncw)]) else numeric(0)
  eta.init <- as.vector(fit$coefficients[(ncw + 1):(ncw + n_eta)])

  # -------------------------
  # 4) baseline hazard increments lambda
  # -------------------------
  M <- nrow(Wtime2)

  Wtime2_phi <- if (ncw > 0) as.vector(Wtime2 %*% phi.init) else rep(0, M)
  linpred <- if (n_eta == 1) as.vector(fixedOrRand.time2) * eta.init else as.vector(fixedOrRand.time2 %*% eta.init)

  temp <- exp(Wtime2_phi + linpred)

  lamb.init <- Index2 / calc_tapply_vect_sum(
    v1 = temp,
    v2 = as.integer(Index1 - 1)
  )

  list(
    phi    = phi.init,
    eta    = eta.init,
    lamb   = lamb.init,
    bBLUP  = bBLUP,
    eta_n  = eta_n,
    B_list = B_list,
    ph     = fit
  )
}
