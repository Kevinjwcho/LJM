############################################################
## Kernel-local posterior LJM: PRES score + Hessian
## (internal helpers used by RefinedfastEM_LLA and RefinedfastEM_LLA_tv)
## - Keep ORIGINAL function names: score_JEL(), PRES_hessian()
## - NO X / NO beta
## - Kernel used ONLY to form local posterior via inv_omega_ker
## - Survival likelihood NOT kernel-weighted (kernel enters only through posterior inputs)
## - var.e score is kernel-weighted (denom is sum of weights, not mi)
##
## Uses your current compiled Rcpp names:
##   calc_mu_surv_noX()
##   calc_Sigma2_surv_noX()
##
## Required in data.mat:
##   Y, Yk, Z, Zk, mi, K, Deltai, W_kernel, w_longK
############################################################


############################################################
## score_JEL()
############################################################
score_JEL <- function(Omega, data.mat, V, b, bhat, Sigmai, S,
                      l0i, l0u, gh.nodes, n, q, nK, nev, Fi, rho = 0){
  
  # ---- Extract fitted values (MLEs) ----
  var.e <- Omega$var.e
  D     <- Omega$D
  eta   <- Omega$eta
  phi <- Omega$phi
  if (is.null(phi)) phi <- numeric(0)
  l0    <- Omega$hazard[, 2]
  ft    <- Omega$hazard[, 1]

  # ---- Extract data objects ----
  Z       <- data.mat$Z
  Zk      <- data.mat$Zk
  Y       <- data.mat$Y
  Yk      <- data.mat$Yk
  mi      <- data.mat$mi

  K       <- data.mat$K
  Deltai  <- data.mat$Deltai
  Di      <- unlist(Deltai)

  # kernel weights
  W_kernel <- data.mat$W_kernel  # stacked weights per subject (len = sum(mi[[i]]))
  w_longK  <- data.mat$w_longK   # list: subject -> marker -> weight vector

  # ---- Gauss-Hermite quadrature ----
  gh <- statmod::gauss.quad.prob(gh.nodes, "normal")
  w  <- gh$weights
  v  <- gh$nodes

  ########################################################
  # Local posterior: kernel only enters here (inv_omega_ker)
  ########################################################
  inv_omega_ker <- calc_inv_omega_kernel(mi, var.e, W_kernel)

  # Mean/Var of (K*phi + eta' b_i | Y_i) using local posterior
  c_list <- data.mat$c
  mu_surv     <- calc_mu_surv_noX(Y, Z, inv_omega_ker, K, D, phi, eta, c_list)
  Sigma2_surv <- calc_Sigma2_surv_noX(Z, inv_omega_ker, D, eta)
  
  # E[exp( linear predictor )] under conditional distribution
  Es_exp <- Esurv_exp(w, v,
                      mu = mu_surv, variance = Sigma2_surv,
                      mu_new = mu_surv, variance_new = Sigma2_surv,
                      l0i = l0i, l0u = l0u, rho = rho)
  
  ########################################################
  # Baseline hazard update pieces (same logic as before)
  ########################################################
  lambda <- update_lambda(Es_exp, l0u, l0, n)
  
  l0.new <- nev / rowSums(do.call(cbind, lambda))
  
  l0u.new <- lapply(l0u, function(x){
    ll <- length(x)
    l0.new[seq_len(ll)]
  })
  
  l0i.new <- numeric(length(Di))
  l0i.new[Di == 0] <- 0
  l0i.new[Di == 1] <- l0.new[match(Fi[Di == 1, 2], ft)]
  
  ########################################################
  # Scores: S(D)
  ########################################################
  Dinv <- solve(D)
  vech.indices <- which(lower.tri(D, diag = TRUE), arr.ind = TRUE)
  dimnames(vech.indices) <- NULL
  
  delta.D <- lapply(seq_len(nrow(vech.indices)), function(d){
    out <- matrix(0, nrow(D), ncol(D))
    ind <- vech.indices[d, 2:1]
    out[ind[1], ind[2]] <- out[ind[2], ind[1]] <- 1
    out
  })
  
  lhs <- sapply(delta.D, function(d) {
    -0.5 * sum(diag(Dinv %*% d))
  })
  
  # score per vech(D) element, aggregated across subjects
  sD <- sapply(seq_along(delta.D), function(iD){
    sum(mapply(function(bi, Si){
      out <- 0.5 * (tcrossprod(bi) + Si) %*% (Dinv %*% delta.D[[iD]] %*% Dinv)
      lhs[iD] + sum(diag(out))
    }, bi = b, Si = Sigmai))
  })
  
  ########################################################
  # Scores: S(var.e)  (kernel-weighted Gaussian score)
  ########################################################
  # Define indices for b blocks (only needed to slice Sigmai into marker blocks)
  if (nK == 1) {
    b.inds <- list(seq_len(q))
  } else {
    b.inds <- split(seq_len(q), cut(seq_len(q), nK, labels = FALSE))
  }
  
  # Build S_list (subject -> marker -> Sigma block) for Sigma.longK
  S_list <- lapply(Sigmai, function(Si) lapply(b.inds, function(ix) Si[ix, ix]))
  
  # Sigma.longK[[i]][[k]] = diag(Zk * S * Zk')
  Sigma.longK <- calc_Sigma_longK(Zk, S_list, nK)
  
  # mu.longK[[i]][[k]] = Yk[,k] - Zk * bhat_block   (no X, no beta)
  mu.longK <- calc_mu_longK_noX(Yk, Zk, bhat, b.inds, nK)
  
  # Ee[[i]][k] = sum_r w_{ikr} ( mu_{ikr}^2 + SigmaDiag_{ikr} )
  Ee <- calc_Ee_kernel(mu.longK, Sigma.longK, w_longK, nK)
  
  # weight sums per subject/marker
  Wsum_mat <- do.call(rbind, lapply(w_longK, function(wi) sapply(wi, sum))) # n x nK
  
  # svar.e summed over subjects: vector length nK
  svar.e_sum <- Reduce("+", lapply(seq_len(n), function(i){
    Ewe_i  <- as.numeric(Ee[[i]])       # length nK
    Wsum_i <- as.numeric(Wsum_mat[i, ]) # length nK
    (-Wsum_i / (2 * var.e)) + (Ewe_i / (2 * (var.e^2)))
  }))
  
  ########################################################
  # Scores: Survival association parameters (eta, phi)
  ########################################################
  Sep_sum <- Setaphi(c(eta, phi),
                     Y, Z, inv_omega_ker,
                     K, D,
                     l0i.new, l0u.new, Di,
                     nK, w, v,
                     mu_surv, Sigma2_surv,
                     eps = 0.0001,
                     c_list = c_list, rho = rho)
  
  ########################################################
  # Full score vector (no beta block)
  ########################################################
  S_out <- c(sD, svar.e_sum, Sep_sum)
  return(S_out)
}



############################################################
## PRES_hessian()
############################################################
PRES_hessian <- function(Omega, data.mat, V, b, bhat, Sigmai, S,
                         l0i, l0u, gh.nodes, n, q, nK, nev, Fi,
                         delta = 0.0001, rho = 0){
  
  # ---- Extract fitted values ----
  D     <- Omega$D
  var.e <- Omega$var.e
  eta   <- Omega$eta
  phi   <- Omega$phi
  
  # ---- Parameter vector (no beta) ----
  params <- c(vech(D), var.e, c(eta, phi))
  len <- length(params)
  
  # ---- Helpers ----
  vech_to_mat <- function(vec, n) {
    mat <- matrix(0, n, n)
    mat[lower.tri(mat, diag = TRUE)] <- vec
    mat <- mat + t(mat)
    diag(mat) <- diag(mat) / 2
    mat
  }
  
  update_Omega <- function(Omega, params) {
    d_dim <- nrow(Omega$D)
    d_len <- d_dim * (d_dim + 1) / 2
    
    # D
    D_new <- vech_to_mat(params[1:d_len], d_dim)
    
    # var.e
    var_e_len   <- length(Omega$var.e)
    var_e_start <- d_len + 1
    var_e_new   <- params[var_e_start:(var_e_start + var_e_len - 1)]
    
    # eta, phi
    eta_len <- length(Omega$eta)
    phi_len <- length(Omega$phi)
    eta_phi_start <- var_e_start + var_e_len
    
    eta_new <- params[eta_phi_start:(eta_phi_start + eta_len - 1)]
    
    if (phi_len > 0) {
      phi_start <- eta_phi_start + eta_len
      phi_new <- params[phi_start:(phi_start + phi_len - 1)]
    } else {
      phi_new <- numeric(0)
    }
    
    Omega$D     <- D_new
    Omega$var.e <- var_e_new
    Omega$eta   <- eta_new
    Omega$phi   <- phi_new
    Omega
  }
  
  # ---- Finite-difference Hessian of score ----
  # Step size (JEL 2.3): `delta` is scaled by the parameter magnitude, and for the
  # entries of D it is halved until all four stencil points keep D positive-definite
  # with a margin (smallest eigenvalue at least half of that of the fitted D). An
  # absolute step of 1e-4 on a near-zero random-slope variance evaluated the score
  # at a singular D and produced Hessian entries of order 1e15, which then ruined the
  # standard errors of eta through the ridged inverse.
  d_dim <- nrow(D); d_len <- d_dim * (d_dim + 1) / 2
  var_e_len <- length(var.e)
  lam_min_D <- min(eigen(D, symmetric = TRUE, only.values = TRUE)$values)
  step_for <- function(i) {
    d_i <- delta * max(1, abs(params[i]))
    if (i <= d_len) {
      pd_ok <- function(d_i) all(vapply(c(-2, -1, 1, 2), function(m) {
        pp <- params; pp[i] <- pp[i] + m * d_i
        ev <- eigen(vech_to_mat(pp[1:d_len], d_dim), symmetric = TRUE, only.values = TRUE)$values
        min(ev) >= 0.5 * lam_min_D
      }, logical(1)))
      for (k in 1:40) { if (pd_ok(d_i)) break; d_i <- d_i / 2 }
    } else if (i <= d_len + var_e_len) {
      d_i <- min(d_i, abs(params[i]) / 4)          # keep sigma^2 - 2*step > 0
    }
    d_i
  }
  H <- matrix(0, nrow = len, ncol = len)
  
  for (i in seq_len(len)) {
    delta_i <- step_for(i)
    params1 <- params2 <- params3 <- params4 <- params
    params1[i] <- params[i] - 2 * delta_i
    params2[i] <- params[i] -     delta_i
    params3[i] <- params[i] +     delta_i
    params4[i] <- params[i] + 2 * delta_i
    
    Omega_new1 <- update_Omega(Omega, params1)
    Omega_new2 <- update_Omega(Omega, params2)
    Omega_new3 <- update_Omega(Omega, params3)
    Omega_new4 <- update_Omega(Omega, params4)
    
    S1 <- score_JEL(Omega_new1, data.mat, V, b, bhat, Sigmai, S,
                    l0i, l0u, gh.nodes, n, q, nK, nev, Fi, rho = rho)
    S2 <- score_JEL(Omega_new2, data.mat, V, b, bhat, Sigmai, S,
                    l0i, l0u, gh.nodes, n, q, nK, nev, Fi, rho = rho)
    S3 <- score_JEL(Omega_new3, data.mat, V, b, bhat, Sigmai, S,
                    l0i, l0u, gh.nodes, n, q, nK, nev, Fi, rho = rho)
    S4 <- score_JEL(Omega_new4, data.mat, V, b, bhat, Sigmai, S,
                    l0i, l0u, gh.nodes, n, q, nK, nev, Fi, rho = rho)

    H[i, ] <- (S1 - 8 * S2 + 8 * S3 - S4) / (12 * delta_i)
  }
  
  H
}