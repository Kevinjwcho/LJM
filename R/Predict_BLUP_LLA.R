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
    
    # Build R_i = diag( sigma_k^2 / w )
    # replicate each sigma_k^2 by number of obs for that k
    m_ik <- sapply(seq_len(nK), function(k) nrow(Zi_blocks[[k]]))
    sig_rep <- rep(var.e, times = m_ik)
    Ri <- diag(sig_rep / (wi + 1e-12))   # avoid division by 0
    
    # closed-form posterior
    A <- Zi %*% D %*% t(Zi) + Ri
    Ainv <- solve(A)
    
    # use (Y - Z c) instead of (Y - X beta)
    resid <- Yi - as.vector(Zi %*% cvec)
    
    mu <- D %*% t(Zi) %*% Ainv %*% resid
    mu_cb <- cvec + mu
    Va <- D - D %*% t(Zi) %*% Ainv %*% Zi %*% D
    
    b_mean[[i]] <- as.numeric(mu_cb)
    b_var[[i]]  <- Va
  }
  
  mean_mat <- do.call(rbind, lapply(b_mean, function(v) matrix(v, nrow=1)))
  colnames(mean_mat) <- paste0(rep(paste0("b",1:nK), each=2), c("_0","_1"))
  
  list(mean = mean_mat, var = b_var)
}