############################################################
## Kernel-local posterior LJM + Time-Varying: PRES score + Hessian
## (internal helper PRES_hessian_LLA_tv used by RefinedfastEM_LLA_tv)
## - NO X / NO beta
## - Kernel used ONLY to form local posterior via inv_omega_ker
## - Time-varying eta via B-spline (etaBBs, BBi, eta.inds)
## - Survival: NOT kernel-weighted (kernel enters only through posterior)
## - var.e score is kernel-weighted
##
## Uses compiled Rcpp names:
##   calc_mu_surv_noX_tv(), calc_Sigma2_surv_noX_tv()
##   etaBBs_compute(), Esurv_exp_t(), update_lambda_tv()
##   Setaphi_lla_tv()
############################################################


############################################################
## score_JEL_LLA_tv()
############################################################
score_JEL_LLA_tv <- function(Omega, data.mat, V, b, bhat, Sigmai, S,
                             l0i, l0u, gh.nodes, n, q, nK, eta.inds, nev, Fi){

  # ---- Extract fitted values (MLEs) ----
  var.e <- Omega$var.e
  D     <- Omega$D
  eta   <- Omega$eta
  phi   <- Omega$phi
  if (is.null(phi)) phi <- numeric(0)
  l0    <- Omega$hazard[, 2]
  ft    <- Omega$hazard[, 1]

  # ---- Safeguard: clamp var.e to positive ----
  var.e <- pmax(var.e, 1e-10)

  # ---- Safeguard: ensure D is positive definite ----
  eig <- eigen(D, symmetric = TRUE)
  if (any(eig$values <= 0)) {
    eig$values <- pmax(eig$values, 1e-10)
    D <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  }

  # ---- Extract data objects ----
  Z       <- data.mat$Z
  Zk      <- data.mat$Zk
  Y       <- data.mat$Y
  Yk      <- data.mat$Yk
  mi      <- data.mat$mi

  K       <- data.mat$K
  Deltai  <- data.mat$Deltai
  Di      <- unlist(Deltai)

  # Kernel weights
  W_kernel <- data.mat$W_kernel
  w_longK  <- data.mat$w_longK

  # B-spline objects
  BBi     <- data.mat$BBi

  # ---- Gauss-Hermite quadrature ----
  gh <- statmod::gauss.quad.prob(gh.nodes, "normal")
  w  <- gh$weights
  v  <- gh$nodes

  ########################################################
  # Local posterior: kernel only enters here (inv_omega_ker)
  ########################################################
  inv_omega_ker <- calc_inv_omega_kernel(mi, var.e, W_kernel)

  # Time-varying etaBBs
  etaBBs <- etaBBs_compute(BBi, eta, eta.inds)

  # Population center c from EM
  c_list <- data.mat$c

  # Mean/Var of survival linear predictor using local posterior (tv)
  mu_surv     <- calc_mu_surv_noX_tv(Y, Z, inv_omega_ker, K, D, phi, etaBBs,
                                     c_list)
  Sigma2_surv <- calc_Sigma2_surv_noX_tv(Z, inv_omega_ker, D, etaBBs)

  # ---- Safeguard: clamp Sigma2_surv to non-negative ----
  Sigma2_surv <- lapply(Sigma2_surv, function(x) pmax(x, 1e-20))

  Es_exp <- Esurv_exp_t(w, v,
                        mu = mu_surv, variance = Sigma2_surv,
                        mu_new = mu_surv, variance_new = Sigma2_surv,
                        l0i = l0i, l0u = l0u)

  ########################################################
  # Baseline hazard update pieces (tv version)
  ########################################################
  lambda <- update_lambda_tv(Es_exp, l0, n)

  l0.new <- nev / rowSums(do.call(cbind, lambda))
  # ---- Safeguard: replace NaN/Inf in l0.new ----
  l0.new[!is.finite(l0.new)] <- 1e-20

  l0u.new <- lapply(l0u, function(x){
    ll <- length(x)
    l0.new[seq_len(ll)]
  })

  l0i.new <- numeric(length(Di))
  l0i.new[Di == 0] <- 0
  l0i.new[Di == 1] <- l0.new[match(Fi[Di == 1, 2], ft)]

  ########################################################
  # Scores: S(D) (same as PRES_MultiJEL.R)
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

  sD <- sapply(seq_along(delta.D), function(iD){
    sum(sapply(b, function(bi){
      out <- 0.5 * tcrossprod(bi) %*% (Dinv %*% delta.D[[iD]] %*% Dinv)
      lhs[iD] + sum(diag(out))
    }))
  })

  ########################################################
  # Scores: S(var.e) (kernel-weighted, same as PRES_MultiJEL.R)
  ########################################################
  if (nK == 1) {
    b.inds <- list(seq_len(q))
  } else {
    b.inds <- split(seq_len(q), cut(seq_len(q), nK, labels = FALSE))
  }

  S_list <- lapply(Sigmai, function(Si) lapply(b.inds, function(ix) Si[ix, ix]))
  Sigma.longK <- calc_Sigma_longK(Zk, S_list, nK)
  mu.longK <- calc_mu_longK_noX(Yk, Zk, bhat, b.inds, nK)
  Ee <- calc_Ee_kernel(mu.longK, Sigma.longK, w_longK, nK)

  Wsum_mat <- do.call(rbind, lapply(w_longK, function(wi) sapply(wi, sum)))

  svar.e_sum <- Reduce("+", lapply(seq_len(n), function(i){
    Ewe_i  <- as.numeric(Ee[[i]])
    Wsum_i <- as.numeric(Wsum_mat[i, ])
    (-Wsum_i / (2 * var.e)) + (Ewe_i / (2 * (var.e^2)))
  }))

  ########################################################
  # Scores: Survival association parameters (eta, phi)
  # Using LLA tv version
  ########################################################
  Sep_sum <- Setaphi_lla_tv(c(eta, phi),
                            Y, Z, inv_omega_ker,
                            K, D,
                            l0i.new, l0u.new, Di,
                            nK, w, v,
                            mu_surv, Sigma2_surv,
                            BBi, eta.inds,
                            eps = 0.0001,
                            c_list = c_list)

  ########################################################
  # Full score vector (no beta block)
  ########################################################
  S_out <- c(sD, svar.e_sum, Sep_sum)
  return(S_out)
}



############################################################
## PRES_hessian_LLA_tv()
############################################################
PRES_hessian_LLA_tv <- function(Omega, data.mat, V, b, bhat, Sigmai, S,
                                l0i, l0u, gh.nodes, n, q, nK, eta.inds, nev, Fi,
                                delta = 0.0001){

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

    D_new <- vech_to_mat(params[1:d_len], d_dim)

    var_e_len   <- length(Omega$var.e)
    var_e_start <- d_len + 1
    var_e_new   <- params[var_e_start:(var_e_start + var_e_len - 1)]

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
  H <- matrix(0, nrow = len, ncol = len)

  for (i in seq_len(len)) {
    params1 <- params2 <- params3 <- params4 <- params
    params1[i] <- params[i] - 2 * delta
    params2[i] <- params[i] -     delta
    params3[i] <- params[i] +     delta
    params4[i] <- params[i] + 2 * delta

    Omega_new1 <- update_Omega(Omega, params1)
    Omega_new2 <- update_Omega(Omega, params2)
    Omega_new3 <- update_Omega(Omega, params3)
    Omega_new4 <- update_Omega(Omega, params4)

    S1 <- score_JEL_LLA_tv(Omega_new1, data.mat, V, b, bhat, Sigmai, S,
                            l0i, l0u, gh.nodes, n, q, nK, eta.inds, nev, Fi)
    S2 <- score_JEL_LLA_tv(Omega_new2, data.mat, V, b, bhat, Sigmai, S,
                            l0i, l0u, gh.nodes, n, q, nK, eta.inds, nev, Fi)
    S3 <- score_JEL_LLA_tv(Omega_new3, data.mat, V, b, bhat, Sigmai, S,
                            l0i, l0u, gh.nodes, n, q, nK, eta.inds, nev, Fi)
    S4 <- score_JEL_LLA_tv(Omega_new4, data.mat, V, b, bhat, Sigmai, S,
                            l0i, l0u, gh.nodes, n, q, nK, eta.inds, nev, Fi)

    H[i, ] <- (S1 - 8 * S2 + 8 * S3 - S4) / (12 * delta)
  }

  H
}
