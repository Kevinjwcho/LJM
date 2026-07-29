
#=============== The DQ Function for Model II ===============#
#=============== Transformation model is fitted for the survival part ===============#

DQfuncGeneric_JEL <- function (model, ptheta, theta, n, Z.st, Y.st, X.st, Ztime, nk, Wtime, Wtime2, Xtime, Xtime2, GQ, Index, Index1, rho, d, wGQ, ncx, ncw, p, ncz, ncz2, b, Ztime2.st, Index0, X, Y, ID, N, Index2, Z, B) { # ptheta means "theta prime"
  ## This might have significant (e-14) numerical difference with the original code. Something odd is happening with VB/muB.
  pbeta <- ptheta$beta
  beta <- theta$beta
  pYsigma2 <- (ptheta$Ysigma) ^ 2
  Ysigma2 <- (theta$Ysigma) ^ 2
  pBsigma <- ptheta$Bsigma
  Bsigma <- theta$Bsigma
  pphi <- ptheta$phi
  phi <- theta$phi
  peta <- ptheta$eta
  eta <- theta$eta
  plamb <- ptheta$lamb
  lamb <- theta$lamb
  n_eta = length(eta)
  
  M <- nrow(Xtime2)
  
  VY <- lapply(1:n, function(i) calc_VY(M = Z.st[[i]], A = Bsigma, b = Ysigma2))  
  VB <-  lapply(1:n, function(i) calc_VB(Bsigma, M2 =  Z.st[[i]], M3 = VY[[i]])) 
  muB <- lapply(1:n, function(i) calc_muB(BSold = Bsigma, Zst = Z.st[[i]], Yst = Y.st[[i]], betaold = beta, VY = VY[[i]], Xst = X.st[[i]]))
  bi.st <- lapply(1:n, function(i) calc_bi_st(v0 = muB[[i]], b, M = VB[[i]]))
  
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
  Wtime2_pphi <- if (ncw > 0) Wtime2 %*% pphi else rep(0, M)
  
  
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
    
    # log.density1 <- log.lamb + as.vector(Wtime_phi.old) + eta.old * Ztime.b # n*GQ matrix #
    log.density1 <- log.lamb + as.vector(Wtime_phi) + eta_Btime.b # n*GQ matrix #
    exp.es <- as.vector(Wtime2_phi) + eta_Btime2.b # M*GQ matrix #  
  } else if (model == "Fixed"){
    # nB <- 0
    log.density1 <- log.lamb + as.vector(Wtime_phi) + eta_time.b # n*GQ matrix #
    exp.es <- as.vector(Wtime2_phi) + eta_time2.b # M*GQ matrix #
  }else{
    stop("Invalid model type")
  } 
  
  const <- matrix(0, n, GQ) # n*GQ matrix #
  
  calc_expM2(exp.es) #A = exp(A)
  
  const[nk != 0, ] <- calc_rowsum_mult((Index), lamb[Index1], exp.es)
  log.density2 <- - log(1 + rho * const) # n*GQ matrix # 
  log.survival <- if (rho > 0) log.density2 / rho else - const # n*GQ matrix # 
  
  f.surv <- exp(d * log.density1 + d * log.density2 + log.survival) # n*GQ matrix #
  deno <- as.vector(f.surv %*% wGQ) # vector of length n #
  Integral <- f.surv / deno # n*GQ matrix f(bi|Oi) #
  CondExp <- (1 + d * rho) / (1 + rho * const) # conditional expectation E(xi|bi,Oi), n*GQ matrix #
  
  # len <- ncx + ncw + p + 2
  len <- ncx + ncw + n_eta + p +  1
  Q <- rep(0, len)
  
  post.bi <- Integral %*% (t(bi) * wGQ) # n*(n*ncz) matrix #
  post.bi <- if(ncz > 1) t(sapply(1:n, function(i) post.bi[i, ((i - 1) * ncz + 1):(i * ncz)])) else 
             matrix(diag(post.bi), nrow = n) # n*ncz matrix #
  
  if (ncz > 1) {    
    tempB <-  fast_rbind_lapply_outerprod( bi.st )  # (n*ncz^2)*GQ matrix #     
  } else {
    tempB <- bi ^ 2
  }
  post.bi2 <- Integral %*% (t(tempB) * wGQ) # n*(n*ncz^2) matrix #
  post.bi2 <- if(ncz > 1) t(sapply(1:n, function(i) post.bi2[i, ((i - 1) * ncz2 + 1):(i * ncz2)])) else 
              matrix(diag(post.bi2), nrow = n) # n*(ncz^2) matrix #
  
  pYmu <- as.vector(X %*% pbeta) + do.call(rbind, lapply(1:n, function(i) Z.st[[i]] %*% bi.st[[i]])) 
  post.resid <- ((Y - pYmu) ^ 2 * Integral[ID, ]) %*% wGQ # vector of length N #
  Q[(ncx + ncw + n_eta + 1)] <- - N / sqrt(pYsigma2) + sum(post.resid) / (pYsigma2 ^ (3 / 2))
  
  pBsigmaInv = solve(pBsigma);
  tempB <- - n * pBsigmaInv / 2 + pBsigmaInv %*% matrix(colSums(post.bi2), ncz, ncz) %*% pBsigmaInv / 2

  ind <- Indexing(ncz)
  Q[(ncx + ncw + n_eta + 2):len] <- as.vector(tapply(c(tempB), ind, sum))

  ## modify end Jun 27
  # There is no Ztime2.b or Ztime.b
  # if (model == 2) {
  #   exp.esp <-as.numeric(Wtime2_pphi) + alpha *Ztime2.b # M*GQ matrix #
  # } else if(model == 1) {
  #   exp.esp <- as.vector(Wtime2_pphi + palpha * Xtime2 %*% pbeta) +  alpha *Ztime2.b # M*GQ matrix #
  #   XZb2 <- as.vector(Xtime2 %*% pbeta) + Ztime2.b # M*GQ matrix #
  # }
  if(model == "TimeVar"){
    exp.esp <-as.numeric(Wtime2_pphi) + eta_Btime2.b # M*GQ matrix #
  }else{
    exp.esp <-as.numeric(Wtime2_pphi) + eta_time2.b # M*GQ matrix #
  }
  
  calc_expM2(exp.esp)

  temp0 <- exp.esp; temp0[1] = temp0[1] + 0 # "touch the variable"
  calc_M1_M2_M3_Hadamard(temp0, CondExp, Integral, as.integer(Index - 1))
  temp1 <- calc_M_v(v = wGQ, M = temp0) #temp1 for lambda I think
  
  #temp0 update to get temp2
  #Since we have serveral etas, we have to make temp0 separately.
  temp0a <- temp0b <- temp0
  if(model == "TimeVar"){
    temp0a <- lapply(1:length(Btime2.b), function(x) temp0*Btime2.b[[x]])
    # temp0 <- Reduce('+', temp0a)
  }else{
    temp0a <- lapply(1:length(time2.b), function(x) temp0*time2.b[[x]])
    # temp0 <- Reduce('+', temp0a)
  }
  
  
  temp2 <- lapply(1:length(temp0a), function(x) calc_M_v(v = wGQ, M = temp0a[[x]])) #compute Mv (matrix vector multiplication)
  post1 <- calc_tapply_vect_sum( temp1, as.integer(Index1 - 1));
  post2 <-lapply(1:length(temp2), function(x) calc_tapply_vect_sum( temp2[[x]], as.integer(Index1 - 1)))  

  # score function for phi
  if (ncw > 0) {
    temp3 <- Wtime2 * temp1 # M*ncw matrix #
    post3 <- as.matrix(apply(temp3, 2, function(x) calc_tapply_vect_sum( x, as.integer(Index1-1))))
    Q[(ncx + 1):(ncx + ncw)] <- colSums(d * Wtime) - colSums(Index2 * post3 / post1) # vector of length ncw #
  }  
  
  # score function for beta and eta previously
  if (model == "TimeVar") {
    sub_Q <- numeric(ncol(B)*nb)
    for(k in 1:ncol(B)){
      for(l in 1:nb){
        sub_Q[nb*(k-1)+l] <- sum(d * post.bi[,l]) - sum(Index2 * post2[[nb*(k-1)+l]] / post1) # for eta
      }
    }
    Q[(ncx + ncw + 1):(ncx + ncw + n_eta)] <- sub_Q
    pResid <- as.vector(Y - X %*% pbeta) - rowSums(Z * post.bi[ID, ]) # vector of length N #
    Q[1:ncx] <- colSums(X * pResid) / pYsigma2 # vector of length ncx # for beta 
  } else { 
    sub_Q <- numeric(nb)
    for(l in 1:nb){
        sub_Q[l] <- sum(d * post.bi[,l]) - sum(Index2 * post2[[l]] / post1) # for eta
    }
    
    Q[(ncx + ncw + 1):(ncx + ncw + n_eta)] <- sub_Q
    pResid <- as.vector(Y - X %*% pbeta) - rowSums(Z * post.bi[ID, ]) # vector of length N #
    Q[1:ncx] <- colSums(X * pResid) / pYsigma2 # vector of length ncx # for beta 
  }
  

return(Q)
}
