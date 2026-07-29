# Standard Error Estimation for Time-Varying JEL (Internal)
#
# This is an internal function for standard error estimation in time-varying JEL.
# @keywords internal
score_JEL_tv <- function(Omega, data.mat, V, b, bhat, Sigmai, S, l0i, l0u, gh.nodes, n, q, nK, eta.inds, nev, Fi){
  # Extract fitted values (MLEs) ----
  beta <- Omega$beta
  var.e <- Omega$var.e
  D <- Omega$D
  eta <- Omega$eta
  phi <- Omega$phi
  l0 <- Omega$hazard[, 2]
  ft <- Omega$hazard[, 1]
  # l0i <- Omega$l0i
  # l0u <- Omega$l0u

  # Extract data objects ----
  # Longitudinal //
  Z <- data.mat$Z; Zk <- data.mat$Zk
  X <- data.mat$X; Xk <- data.mat$Xk
  Y <- data.mat$Y; Yk <- data.mat$Yk
  mi <- data.mat$mi

  # Survival //
  K <- data.mat$K
  Deltai <- data.mat$Deltai
  Di = unlist(Deltai)
  BBi <- data.mat$BBi

  gh <- statmod::gauss.quad.prob(gh.nodes, 'normal')
  w <- gh$weights; v <- gh$nodes

  # Hazard ----
  inv_omega <- calc_inv_omega(mi, var.e)
  etaBBs <- etaBBs_compute(BBi, eta, eta.inds)
  mu_surv <- calc_mu_surv_tv(X, Y, Z, inv_omega, K, D, beta, phi, etaBBs, nK)
  Sigma2_surv <- calc_Sigma2_surv_tv(X, Z, inv_omega, D, etaBBs)

  Es_exp <- Esurv_exp_t(w, v, mu = mu_surv, variance = Sigma2_surv,
                        mu_new = mu_surv, variance_new = Sigma2_surv, l0i, l0u)

  # The baseline hazard, \lambda ----
  lambda <- update_lambda_tv(Es_exp, l0, n)

  l0.new <- nev/rowSums(do.call(cbind, lambda))
  l0u.new <- lapply(l0u, function(x){
    ll <- length(x); l0.new[1:ll]
  })
  l0i.new <- c()
  l0i.new[which(Di == 0)] <- 0
  l0i.new[which(Di == 1)] <- l0.new[match(Fi[which(Di==1), 2], ft)]

  # (\eta, \phi) ----

  # Scores ----------------------------------------------------------------
  # S(D) ----

  Dinv <- solve(D)
  vech.indices <- which(lower.tri(D, diag = T), arr.ind = T)
  dimnames(vech.indices) <- NULL
  delta.D <- lapply(1:nrow(vech.indices), function(d){
    out <- matrix(0, nrow(D), ncol(D))
    ind <- vech.indices[d, 2:1]
    out[ind[1], ind[2]] <- out[ind[2], ind[1]] <- 1 # dD/dvech(d)_i
    out
  })

  lhs <- sapply(delta.D, function(d) {
    -0.5 * sum(diag(Dinv %*% d))
  })

  sDi <- function(i) {
    mapply(function(b) {
      out <- 0.5 * tcrossprod(b) %*% (Dinv %*% delta.D[[i]] %*% Dinv)   # (Sigmai + tcrossprod(b))?
      lhs[i] + sum(diag(out))
    },
    b = b,
    SIMPLIFY = T)
  }

  sD <- sapply(1:nrow(vech.indices), sDi)
  sD <- lapply(1:nrow(sD), function(x) sD[x, ]) # Cast to list
  sD_sum = Reduce('+', sD)

  # Longitudinal ----
  # Define indices for beta and b
  if(nK == 1){
    b.inds <- list(seq(q))
    beta.inds <- list(seq(length(beta)))
  }else{
    b.inds <- split(seq(q), cut(seq_along(seq(q)), nK, labels = F))
    beta.inds <- split(seq(length(beta)), cut(seq_along(seq(length(beta))), nK, labels = F))
  }


  # S(\beta) -----

  sbeta <- mapply(function(X, V, Y, Z, b){
    crossprod(X, solve(V) %*% (Y - X %*% beta - Z %*% b))
  },
  X = X, V = V, Y = Y, Z = Z, b = b, SIMPLIFY = F)
  sbeta_sum = Reduce('+', sbeta)

  # S(var.e) ----
  svar.e <- list()
  for(i in 1:n){
    temp <- numeric(nK)
    for(k in 1:nK){
      tr_sigma_tilde <- sum(diag(Zk[[i]][[k]] %*% Sigmai[[i]][b.inds[[k]], b.inds[[k]]] %*% t(Zk[[i]][[k]])))
      mu_tilde <- crossprod(Yk[[i]][, k] - Xk[[i]][[k]] %*% beta[beta.inds[[k]]] -
                              Zk[[i]][[k]] %*% b[[i]][b.inds[[k]]])

      temp[k] <- -mi[[i]][k]/(2 * var.e[k]) + 1/(2 * var.e[k]^2) * (mu_tilde  + tr_sigma_tilde)
    }
    svar.e[[i]] <- temp
  }
  svar.e_sum = Reduce('+', svar.e)

  # Survival parameters ----
  # inv_omega <- calc_inv_omega(mi, var.e)
  # etaBBs <- etaBBs_compute(BBi, eta, eta.inds)
  # mu_surv <- calc_mu_surv(X, Y, Z, inv_omega, K, D, beta, phi, etaBBs, nK)
  # Sigma2_surv <- calc_Sigma2_surv(X, Z, inv_omega, D, etaBBs)

  Sep_sum = Setaphi_t(c(eta, phi), X, Y, Z, inv_omega, K, D, beta,
                    l0i.new, l0u.new, Di, nK, w, v, mu_surv, Sigma2_surv, BBi, eta.inds, eps = 0.0001)

  # Forming information -----------------------------------------------------

  S <- c(sD_sum, sbeta_sum, svar.e_sum, Sep_sum)
  # NB RHS should = 0, however due to approximate nature of the approach (i.e. Omega aren't MLEs), we leave this term in.
  return(S)
}



# Hessian -----------------------------------------------------------------

PRES_hessian_tv = function(Omega, data.mat, V, b, bhat, Sigmai, S, l0i, l0u, gh.nodes, n, q, nK, eta.inds, nev, Fi, delta = 0.0001){
  # Extract fitted values (MLEs) ----
  beta <- Omega$beta
  var.e <- Omega$var.e
  D <- Omega$D
  eta <- Omega$eta
  phi <- Omega$phi

  # parameters
  params = c(vech(D), beta, var.e, c(eta, phi))
  len = length(params)

  vech_to_mat <- function(vec, n) {
    mat <- matrix(0, n, n)
    mat[lower.tri(mat, diag = TRUE)] <- vec
    mat <- mat + t(mat)
    diag(mat) <- diag(mat) / 2
    return(mat)
  }

  update_Omega <- function(Omega, params) {
    # Determine size of D
    d_dim <- nrow(Omega$D)

    # Extract elements for D based on the triangle
    d_len = (d_dim * (d_dim + 1) / 2)
    D_new <- vech_to_mat(params[1:(d_dim * (d_dim + 1) / 2)], d_dim)

    # Extract beta
    beta_len <- length(Omega$beta)
    beta_start <- d_len + 1
    beta_new <- params[beta_start:(beta_start + beta_len - 1)]

    # Extract var.e
    var_e_len <- length(Omega$var.e)
    var_e_start <- beta_start + beta_len
    var_e_new <- params[var_e_start:(var_e_start + var_e_len - 1)]

    # Extract eta and phi
    eta_len <- length(Omega$eta)
    phi_len <- length(Omega$phi)
    eta_phi_start <- beta_start + beta_len + var_e_len
    eta_new <- params[eta_phi_start:(eta_phi_start + eta_len - 1)]
    phi_new <- params[(eta_phi_start + eta_len):(eta_phi_start + eta_len + phi_len - 1)]

    # Update Omega
    Omega$D <- D_new
    Omega$beta <- beta_new
    Omega$var.e <- var_e_new
    Omega$eta <- eta_new
    Omega$phi <- phi_new

    return(Omega)
  }


  H = matrix(0, nrow = len, ncol = len)
  for(i in 1:len){
    params1 <- params2 <- params3 <- params4 <- params
    params1[i] = params[i] - 2* delta
    params2[i] = params[i] - delta
    params3[i] = params[i] + delta
    params4[i] = params[i] + 2*delta
    Omega_new1 = update_Omega(Omega, params1)
    Omega_new2 = update_Omega(Omega, params2)
    Omega_new3 = update_Omega(Omega, params3)
    Omega_new4 = update_Omega(Omega, params4)
    S1 = score_JEL_tv(Omega_new1, data.mat, V, b, bhat, Sigmai, S, l0i, l0u, gh.nodes, n, q, nK, eta.inds, nev, Fi)
    S2 = score_JEL_tv(Omega_new2, data.mat, V, b, bhat, Sigmai, S, l0i, l0u, gh.nodes, n, q, nK, eta.inds, nev, Fi)
    S3 = score_JEL_tv(Omega_new3, data.mat, V, b, bhat, Sigmai, S, l0i, l0u, gh.nodes, n, q, nK, eta.inds, nev, Fi)
    S4 = score_JEL_tv(Omega_new4, data.mat, V, b, bhat, Sigmai, S, l0i, l0u, gh.nodes, n, q, nK, eta.inds, nev, Fi)

    H[i, ] <- (S1 - 8 * S2 + 8 * S3 - S4) / (12 * delta)
  }

  # names(H) <- names(params)
  return(H)
}
