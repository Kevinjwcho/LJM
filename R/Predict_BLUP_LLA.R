Predict_BLUP_LLA <- function(prep_te, fitLLAJEL){
  
  D     <- fitLLAJEL$coefficients$D
  var.e <- fitLLAJEL$coefficients$var.e
  nK    <- length(var.e)
  
  # subject count (assumes prep_te$ID[[k]] indexes 1..nsub consistently)
  nsub <- length(prep_te$ni[[1]])
  
  # build stacked Z_i, Y_i, w_i for each subject
  Z.st <- lapply(seq_len(nK), function(k){
    split(prep_te$Z[[k]], prep_te$ID[[k]]) |>
      lapply(function(v) matrix(v, ncol = ncol(prep_te$Z[[k]])))
  })
  Y.st <- lapply(seq_len(nK), function(k){
    split(prep_te$Y[[k]], prep_te$ID[[k]])
  })
  W.st <- lapply(seq_len(nK), function(k){
    split(prep_te$W_ker[[k]], prep_te$ID[[k]])
  })
  
  # centering c (length 2*nK): use the fitted center of the random-effects
  # distribution (paper BLUP uses the model parameter c). Fall back to the
  # prep-side kernel-LLA center for backward compatibility if c is absent.
  c_fit <- fitLLAJEL$coefficients$c
  if (!is.null(c_fit)) {
    cvec <- as.numeric(c_fit)
  } else {
    cvec <- as.numeric(unlist(prep_te$cLLA))
  }
  
  b_mean <- vector("list", nsub)
  b_var  <- vector("list", nsub)
  
  for (i in seq_len(nsub)){
    
    Zi_blocks <- lapply(seq_len(nK), function(k) Z.st[[k]][[i]])
    Zi <- as.matrix(Matrix::bdiag(Zi_blocks))
    
    Yi <- unlist(lapply(seq_len(nK), function(k) Y.st[[k]][[i]]))
    wi <- unlist(lapply(seq_len(nK), function(k) W.st[[k]][[i]]))
    
    # Precision (information) form of the kernel-weighted BLUP / MAP predictor:
    #   b_hat = (D^-1 + Z' Omega Z)^-1 (D^-1 c + Z' Omega y),   Omega = diag(w / sigma_k^2)
    # This is the manuscript's Section 3.1 expression.  It is algebraically
    # identical to the covariance form  c + D Z'(Omega^-1 + Z D Z')^-1 (y - Z c)
    # used previously (Woodbury), but inverts a 2K x 2K matrix instead of an
    # N_i x N_i one and never forms Omega^-1, so zero-weight observations no
    # longer enter as sigma^2/1e-12 (~1e12) diagonal entries.
    m_ik <- sapply(seq_len(nK), function(k) nrow(Zi_blocks[[k]]))
    sig_rep <- rep(var.e, times = m_ik)
    om <- wi / sig_rep                    # Omega_i diagonal; zero weight -> zero information
    ZtO <- crossprod(Zi, Zi * om)         # Z' Omega Z   (2K x 2K)
    Dinv <- solve(D)
    Ai <- Dinv + ZtO
    Ai <- 0.5 * (Ai + t(Ai))              # enforce symmetry against round-off
    rhs <- as.vector(Dinv %*% cvec) + as.vector(crossprod(Zi, Yi * om))
    mu_cb <- as.numeric(solve(Ai, rhs))
    Va <- solve(Ai)
    
    b_mean[[i]] <- as.numeric(mu_cb)
    b_var[[i]]  <- Va
  }
  
  mean_mat <- do.call(rbind, lapply(b_mean, function(v) matrix(v, nrow=1)))
  colnames(mean_mat) <- paste0(rep(paste0("b",1:nK), each=2), c("_0","_1"))
  
  list(mean = mean_mat, var = b_var)
}