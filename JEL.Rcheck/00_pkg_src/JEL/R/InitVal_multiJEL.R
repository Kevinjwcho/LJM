
#=============== Initial Value Calculation for JEL ===============#

InitVal_multiJEL <- function (beta, model, n, n_LME, X, Z, Y, bBLUP, ID, Xtime, Ztime, Xtime2, Ztime2, Indcs, start, event, stop, W , ncw, Wtime2, nk, Wtime22, d, Wtime, ID_surv){
  
  Index = Indcs$Index  
  Index0 = Indcs$Index0
  Index1 = Indcs$Index1
  Index2 = Indcs$Index2

  
  M <- lapply(Xtime2, nrow)
  nrow_surv <- nrow(W)
  
  n_eta <- sapply(bBLUP, ncol) %>% sum
  fixedOrRand.time <- lapply(1:n_LME, function(x) bBLUP[[x]][ID_surv,]) # vector of length of lmm model #
  fixedOrRand.time <- do.call('cbind', fixedOrRand.time)
  fixedOrRand.time2 <- lapply(1:n_LME, function(x) bBLUP[[x]][Index,]) # vector of length of lmm model #
  fixedOrRand.time2 <- do.call('cbind', fixedOrRand.time2)
  
  eta_n <- c()
  for(i in 1:n_LME){
    if(ncol(bBLUP[[i]]) == 1){
      eta_ni <- paste0("eta_b", i, "_0")
    }else{
      eta_ni <- paste0("eta_b", i, "_", c(0:(ncol(bBLUP[[i]])-1)))  %>% as.vector
    }
    eta_n <- c(eta_n, eta_ni)
  }
  
  
  
  
  #========== first fit the Cox model ==========#
  data.init <- data.frame(start = start, stop = stop, event = event, W = W, fixedOrRand = fixedOrRand.time)
  fit <- if (ncw > 0) coxph(Surv(start, stop, event) ~ W + fixedOrRand.time, data = data.init) else coxph(Surv(start, stop, event) ~ fixedOrRand.time, data = data.init)
  phi.old <- if (ncw > 0) fit$coefficients[1:ncw] else numeric(0)
  
  eta.old <- fit$coefficients[(ncw+1):(ncw+n_eta)] %>% as.vector
  
  
  
  Wtime2_phi.old <- if (ncw > 0) Wtime2 %*% phi.old else rep(0, M)
  if(n_eta == 1){
    temp <- as.vector(exp(Wtime2_phi.old + fixedOrRand.time2*eta.old)) # M*1 vector # 
  } else{
    temp <- as.vector(exp(Wtime2_phi.old + fixedOrRand.time2%*%eta.old)) # M*1 vector #  
  }
  
  lamb.old <- Index2 / calc_tapply_vect_sum(  v1 = temp, v2 = as.integer(Index1 - 1)) # vector of length n_u # Index 1 with length M.
  
  
  phi.new <- phi.old
  eta.new <- eta.old
  lamb.new <- lamb.old
  
  B <- NULL
  result <- list(phi = phi.new, eta = eta.new, lamb = lamb.new, B=B, eta_n = eta_n, ph = fit)
  return(result)
}
