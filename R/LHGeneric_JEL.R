
#=============== Function to Calculate the Likelihood Value for Model II ===============#
#LHGeneric <- function (theta){
LHGeneric_JEL <- function (theta, n, Z.st, Y.st, X.st, b, Ztime, nk, Index0, model, Wtime, Xtime, Wtime2, Xtime2, GQ, Index, Index1, rho, d, wGQ, ncz, Ztime2.st, ncw, B) {
  
  beta <- theta$beta
  Ysigma2 <- (theta$Ysigma) ^ 2
  Bsigma <- theta$Bsigma
  phi <- theta$phi
  eta <- theta$eta
  lamb <- theta$lamb
  
  M <- nrow(Xtime2)
  
  VY <- lapply(1 : n, function(i) calc_VY(M = Z.st[[i]], A = Bsigma, b = Ysigma2))  
  VB <-  lapply(1 : n, function(i) calc_VB(Bsigma, M2 =  Z.st[[i]], M3 = VY[[i]])) 
  muB <- lapply(1 : n, function(i) calc_muB(BSold = Bsigma, Zst = Z.st[[i]], Yst = Y.st[[i]], betaold = beta, VY = VY[[i]], Xst = X.st[[i]]))
  bi.st <- lapply(1 : n, function(i) calc_bi_st(v0 = muB[[i]], b, M = VB[[i]]))
  
  bi <- do.call(rbind, bi.st)
  
  # Ztime.b <- do.call(rbind, lapply(1 : n, function(i) Ztime[i, ] %*% bi.st[[i]])) # n*GQ matrix #
  # Ztime2.b <-fast_lapply_length(Ztime2.st, bi.st, (1 : n)[nk != 0] - 1)# M*GQ matrix #
  
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
  
  
  log.lamb <- log(lamb[Index0])
  log.lamb[is.na(log.lamb)] <- 0
  
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
    
    # log.density1 <- log.lamb + as.vector(Wtime_phi) + eta * Ztime.b # n*GQ matrix #
    log.density1 <- log.lamb + as.vector(Wtime_phi) + eta_Btime.b # n*GQ matrix #
    h.s <- as.vector(Wtime2_phi) + eta_Btime2.b # M*GQ matrix #  
  } else if (model == "Fixed"){
    # nB <- 0
    log.density1 <- log.lamb + as.vector(Wtime_phi) + eta_time.b # n*GQ matrix #
    h.s <- as.vector(Wtime2_phi) + eta_time2.b # M*GQ matrix #
  }else{
    stop("Invalid model type")
  } 
  
  const <- matrix(0, n, GQ) # n*GQ matrix #
  
  calc_expM2(h.s) #A = exp(A)
  temp0a <- h.s * lamb[Index1]
  
  const[nk != 0, ] <- calc_rowsum((Index), temp0a)
  log.density2 <- - log(1 + rho * const) # n*GQ matrix # 
  log.survival <- if (rho > 0) log.density2 / rho else - const # n*GQ matrix # 
  
  f.surv <- exp(d * log.density1 + d * log.density2 + log.survival) # n*GQ matrix #
  deno <- as.vector(f.surv %*% wGQ) # vector of length n #
  Integral <- f.surv / deno # n*GQ matrix f(bi|Oi) #
  
  f.long <- sapply(1 : n, function(i) calc_MVND(Y.st[[i]], as.vector(X.st[[i]] %*% beta), VY[[i]]))
  
  # lgLik <- sum(log(f.long * deno / (pi ^ (ncz / 2)))) # pi^ related to sigma^2
  # vector of length n #
  return(sum(log(f.long * deno / (pi ^ (ncz / 2)))))
}
