Predict_BLUP_Multi <- function(fit.JEL, fit.LMM_list, LME_dat.te){
  
  #### For test
  # fit.JEL = JEL_models[[1]]$JELmodels[[1]]
  # fit.LMM_list = JEL_models[[1]]$LMMmodels
  
  
  ### real start
  sigma_e <- fit.JEL$coefficients$Ysigma
  Bsigma <- fit.JEL$coefficients$Bsigma
  beta_est <- fit.JEL$coefficients$beta
  
  
  # formLongZ <- formula(fit.LMM$modelStruct$reStruct[[1]]) # random effect formula
  formLongZ <- lapply(fit.LMM_list, function(x) formula(x$modelStruct$reStruct[[1]]))
  # mfLongZ <- model.frame(terms(formLongZ), data = LME_dat.te) # Z
  mfLongZ <- lapply(formLongZ, function(x) model.frame(terms(x), data = LME_dat.te))
  # TermsLongZ <- attr(mfLongZ, "terms") 
  TermsLongZ <- lapply(mfLongZ, function(x) attr(x, "terms"))
  # Z <- as.data.frame(model.matrix(formLongZ, mfLongZ)) 
  Z <- lapply(1:length(fit.LMM_list), function(i) as.data.frame(model.matrix(formLongZ[[i]], mfLongZ[[i]])))
  # split the matrix by id
  Z_new <- lapply(Z, function(x) split(x, LME_dat.te$id))
  Z_new <- lapply(Z_new, function(x) lapply(x, as.matrix))
  n_te <- lapply(Z_new, function(x) length(x))
  if(length(unique(n_te)) > 1){
    stop("The number of observations for each model is not the same")
  }
  Z_bd = lapply(1:n_te[[1]], function(name) {
    sub = sapply(Z_new, function(sublist) sublist[[name]], simplify = FALSE)
    return(as.matrix(Matrix::bdiag(sub)))
  })
  
  
  
  TermsLongQ <- lapply(fit.LMM_list, function(x) x$terms)
  formLongQ <- lapply(fit.LMM_list, function(x) formula(x))
  mfLongQ <- lapply(1:length(fit.LMM_list), function(i) model.frame(TermsLongQ[[i]], data = LME_dat.te))
  
  Q <- lapply(1:length(fit.LMM_list), function(i) as.data.frame(model.matrix(formLongQ[[i]], mfLongQ[[i]])))
  Q_new <- lapply(Q, function(x) split(x, LME_dat.te$id))
  Q_new <- lapply(Q_new, function(x) lapply(x, as.matrix))
  Q_bd = lapply(1:n_te[[1]], function(name) {
    sub = sapply(Q_new, function(sublist) sublist[[name]], simplify = FALSE)
    return(as.matrix(Matrix::bdiag(sub)))
  })
  
  Y <- lapply(1:length(fit.LMM_list), function(i) as.vector(model.response(mfLongQ[[i]], "numeric")))
  # Y_new <- split(Y, LME_dat.te$id)
  Y_new <- lapply(Y, function(x) split(x, LME_dat.te$id))
  Y_vec = lapply(1:n_te[[1]], function(name) {
    sub = sapply(Y_new, function(sublist) sublist[[name]], simplify = FALSE)
    return(as.vector(unlist(sub)))
  })
  
  # split Bsigma into the square matrix with the dimension lapply(Z, ncol)
  # Bsigma_list = list()
  # beta_list = list()
  # n_b = sapply(Z, ncol)
  # for(i in 1:length(fit.LMM_list)){
  #   nb = n_b[i]
  #   Bsigma_list[[i]] = matrix(Bsigma[(1:(nb)) + (i-1)*nb], nb, nb)
  #   beta_list[[i]] = beta_est[(1:(nb)) + (i-1)*nb]
  # }
  
  nk = lapply(Z_new, function(x) sapply(x, nrow))
  nk_vec = lapply(1:n_te[[1]], function(name) {
    sub = sapply(nk, function(sublist) sublist[[name]], simplify = FALSE)
    return(as.vector(unlist(sub)))
  })
  sigma.e.list = lapply(nk_vec, function(x) diag(rep(sigma_e, x)))
  
  BLUP_new <- list()
  if(is.matrix(Bsigma)){
    # j is the number of observations for the i-th model
    BLUP_new = list(mean = sapply(1:n_te[[1]], function(j) ((Bsigma%*%t(Z_bd[[j]]))%*%solve(sigma.e.list[[j]]*diag(rep(1, nrow(Z_bd[[j]]))) + Z_bd[[j]]%*%Bsigma %*% t(Z_bd[[j]]))%*%(Y_vec[[j]] - Q_bd[[j]]%*%beta_est))),
                    var = lapply(1:n_te[[1]], function(j) Bsigma - Bsigma%*%t(Z_bd[[j]])%*%solve(sigma.e.list[[j]]*diag(rep(1, nrow(Z_bd[[j]]))) + Z_bd[[j]]%*%Bsigma %*% t(Z_bd[[j]]))%*%Z_bd[[j]]%*%t(Bsigma))
    )
  }else{
    stop("Bsigma is not a matrix")
  }
  
  # BLUP_new <- list("mean" = EBLUP_new, "var" = varBLUP_new)
  return(BLUP_new)
  
}
