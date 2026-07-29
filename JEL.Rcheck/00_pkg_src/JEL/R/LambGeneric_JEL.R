
#========== Function to Obtain Lamb Given Other Finite Dimensional Parameters for Model II ==========#
#=============== Transformation model is fitted for the survival part ===============#

LambGeneric_JEL <- function (para, lamb.init, tol, iter, ncz, ncx, ncw, n_eta, n, Z.st, Y.st, X.st, b, Ztime, Ztime2.st, nk, Wtime, Xtime, Wtime2, Xtime2, rho, Index0, Index1, Index, wGQ, model, GQ, d, Index2, B){

  para.list <- Vec2List_JEL(para, ncx, ncz, ncw, n_eta)
  beta <- para.list$beta
  Ysigma2 <- (para.list$Ysigma) ^ 2
  Bsigma <- para.list$Bsigma
  phi <- para.list$phi
  eta <- para.list$eta
  
  
  M <- nrow(Xtime2)
  
  VY <- lapply(1 : n, function(i) calc_VY(M = Z.st[[i]], A = Bsigma, b = Ysigma2))  
  VB <-  lapply(1 : n, function(i) calc_VB(Bsigma, M2 =  Z.st[[i]], M3 = VY[[i]])) 
  muB <- lapply(1 : n, function(i) calc_muB(BSold = Bsigma, Zst = Z.st[[i]], Yst = Y.st[[i]], betaold = beta, VY = VY[[i]], Xst = X.st[[i]]))
  bi.st <- lapply(1 : n, function(i) calc_bi_st(v0 = muB[[i]], b, M = VB[[i]]))
  
  bi <- do.call(rbind, bi.st)
  
  nb <- nrow(bi.st[[1]])
  
  if(model == "Fixed"){
    time.b = list()
    for(i in 1:nb){
      time.b[[i]] <- do.call(rbind, lapply(1:n, function(j) bi.st[[j]][i,]))
    }
    time2.b <- lapply(1:nb, function(i) time.b[[i]][Index, ])
    temp.eta_time.b <- lapply(1:nb, function(i) time.b[[i]]*eta[i])
    eta_time.b <- Reduce('+', temp.eta_time.b)
    temp.eta_time2.b <- lapply(1:nb, function(i) time2.b[[i]]*eta[i])
    eta_time2.b <- Reduce('+', temp.eta_time2.b)
  }
  
  
  # log.lamb <- log(lamb[Index0])
  # log.lamb[is.na(log.lamb)] <- 0
  
  Wtime_phi <- if (ncw > 0) Wtime %*% phi else rep(0, n)
  Wtime2_phi <- if (ncw > 0) Wtime2 %*% phi else rep(0, M)
  
  
  ## log likelihood 
  if (model == "TimeVar"){
    # nb <- nrow(bi.st[[1]])
    nB <- ncol(B)
    
    b_list <- list()
    b_list2 <- list()
    for(j in 1:nb){
      b_list[[j]] <- sapply(1:n, function(i) bi.st[[i]][j, ]) %>% t()
      b_list2[[j]] <- sapply(Index, function(i) bi.st[[i]][j, ]) %>% t()  
    }
    
    
    Btime.b <- list()
    for(i in 1:nB){
      for(j in 1:nb){
        list_ind <- nb*(i-1)+j
        B_diag <- diag(B[,i])
        Btime.b[[list_ind]] <- B_diag%*%b_list[[j]]
      }
    }
    
    B2 <- B[Index, ]
    if(nB == 1) B2 <- matrix(B2, ncol = nB)
    
    Btime2.b <- list()
    for(i in 1:nB){
      for(j in 1:nb){
        list_ind <- nb*(i-1)+j
        B2_diag <- diag(B2[,i])
        Btime2.b[[list_ind]] <- B2_diag%*%b_list2[[j]]
      }
    }
    
    
    eta_Btime.b_list <- lapply(1:length(Btime.b), function(i) Btime.b[[i]] * eta[i])
    eta_Btime.b <- Reduce('+', eta_Btime.b_list)
    
    eta_Btime2.b_list <- lapply(1:length(Btime2.b), function(i) Btime2.b[[i]] * eta[i])
    eta_Btime2.b <- Reduce('+', eta_Btime2.b_list)
    
    
    eta.h <- as.vector(Wtime_phi) + eta_Btime.b 
    exp.es <- as.vector(Wtime2_phi) + eta_Btime2.b # M*GQ matrix #  
  } else if (model == "Fixed"){
    eta.h <- as.vector(Wtime_phi) + eta_time.b 
    exp.es <- as.vector(Wtime2_phi) + eta_time2.b # M*GQ matrix #  
  }else{
    stop("Invalid model type")
  } 
  
  # 
  # M <- nrow(Xtime2)
  # 
  # VY <- lapply(1:n, function(i) calc_VY(M = Z.st[[i]], A = Bsigma, b = Ysigma2 ))  
  # VB <-  lapply(1:n, function(i) calc_VB( Bsigma ,M2 =  Z.st[[i]], M3 = VY[[i]])) 
  # muB <- lapply(1:n, function(i) calc_muB( BSold=Bsigma , Zst=Z.st[[i]], Yst=Y.st[[i]], betaold=beta ,VY= VY[[i]], Xst=X.st[[i]]))
  # bi.st <- lapply(1:n, function(i) calc_bi_st(v0=muB[[i]], b ,M = VB[[i]]) ) 
  # 
  # # each element is ncz*GQ matrix #
  # bi <- do.call(rbind, bi.st) # (n*ncz)*GQ mat rix #
  # Ztime.b <- do.call(rbind, lapply(1:n, function(i) Ztime[i, ] %*% bi.st[[i]])) # n*GQ matrix #
  # Ztime2.b <-fast_lapply_length(Ztime2.st, bi.st, (1:n)[nk != 0] - 1) # M*GQ matrix #  
  # 
  # Wtime_phi <- if (ncw > 0) Wtime %*% phi else rep(0, n)
  # Wtime2_phi <- if (ncw > 0) Wtime2 %*% phi else rep(0, M)
  # 
  # if (model == 2) { 
  #   eta.h <- as.vector(Wtime_phi) + alpha * Ztime.b # n*GQ matrix #
  # } else if(model ==1) {
  #   eta.h <- as.vector(Wtime_phi + alpha * Xtime %*% beta) + alpha * Ztime.b # n*GQ matrix #
  # } else {
  #   stop("Invalid model type")
  # }
  # calc_v_a( Ztime2.b, alpha); # Ztime2.b gets altered
  # if ( model == 2) {
  #   exp.es<- as.numeric(Wtime2_phi) + Ztime2.b  
  # } else {
  #   exp.es<- as.numeric(Wtime2_phi + alpha * Xtime2 %*% beta) + Ztime2.b  
  # }
  calc_expM2(exp.es)
  lamb.old <- lamb.init
  err <- 1
  
  for (step in 1:iter) {
    Lamb.old <- cumsum(lamb.old)
    
    log.lamb <- log(lamb.old[Index0])
    log.lamb[is.na(log.lamb)] <- 0
    log.density1 <- log.lamb + eta.h # n*GQ matrix #
    const <- matrix(0, n, GQ) # n*GQ matrix 
    const[nk != 0, ] <- calc_rowsum_mult((Index), lamb.old[Index1], exp.es)
    log.density2 <- - log(1 + rho * const) # n*GQ matrix # 
    log.survival <- if(rho > 0) log.density2 / rho else - const # n*GQ matrix # 
    
    f.surv <- exp(d * log.density1 + d * log.density2 + log.survival) # n*GQ matrix #
    deno <- as.vector(f.surv %*% wGQ) # vector of length n #
    Integral <- f.surv / deno # n*GQ matrix #
    CondExp <- (1 + d * rho) / (1 + rho * const) # conditional expectation E(xi|bi,Oi), n*GQ matrix #
    
    # WE DO IN PLACE MULTIPLICATION / variable 'tempLamb0' is holding the results    
    tempLamb0 <- exp.es; tempLamb0[1] = tempLamb0[1] +0 # "touch the variable"
    calc_M1_M2_M3_Hadamard(tempLamb0, CondExp ,  Integral, as.integer(Index-1))
    tempLamb <- calc_M_v(v =wGQ, M =tempLamb0)
    postLamb <- calc_tapply_vect_sum( v1= tempLamb, v2= as.integer(Index1-1)); ## Check this!
    lamb.new <- Index2 / postLamb[postLamb!=0]
    
    Lamb.new <- cumsum(lamb.new)
    err <- max(abs((Lamb.new - Lamb.old) / Lamb.old))
    if (step > 3 & err < tol) break
    lamb.old <- lamb.new
  }
  converge <- as.numeric(step <= iter & err < tol)
  return(list(lamb = lamb.new, converge = converge))
}
