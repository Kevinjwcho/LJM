
#=============== EM iteration Using Adaptive Gaussian Quadrature for Model I&II ===============#
#=============== Transformation model is fitted for the survival part ===============#

EMiterJEL <- function (theta.old, n, Z.st, Ztime, Ztime2.st, nk, Wtime2, Xtime2, GQ, rho, wGQ, d, Y.st, X.st, ncz, ncz2, b, model, Wtime, Xtime, X, Y, ID, N, ncw, Wtime22, ncx, Xtime22, Z, X2.sum, Indcs, B, cntrlLst){
  # Use apply instead of matrix calculation #
  
  # Get Old Estimates #
  beta.old <- theta.old$beta
  Ysigma2.old <- (theta.old$Ysigma) ^ 2
  Bsigma.old <- theta.old$Bsigma
  phi.old <- theta.old$phi
  eta.old <- theta.old$eta
  
  Index = Indcs$Index  
  Index0 = Indcs$Index0
  Index1 = Indcs$Index1
  Index2 = Indcs$Index2
  
  lamb.old <- theta.old$lamb
  
  M <- nrow(Xtime2)
  
  VY <- lapply(1 : n, function(i) calc_VY(M = Z.st[[i]], A = Bsigma.old, b = Ysigma2.old))  
  VB <-  lapply(1 : n, function(i) calc_VB(Bsigma.old, M2 =  Z.st[[i]], M3 = VY[[i]])) 
  muB <- lapply(1 : n, function(i) calc_muB(BSold = Bsigma.old, Zst = Z.st[[i]], Yst = Y.st[[i]], betaold = beta.old, VY = VY[[i]], Xst = X.st[[i]]))
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
    temp.eta_time.b <- lapply(1:nb, function(i) time.b[[i]]*eta.old[i])
    eta_time.b <- Reduce('+', temp.eta_time.b)
    temp.eta_time2.b <- lapply(1:nb, function(i) time2.b[[i]]*eta.old[i])
    eta_time2.b <- Reduce('+', temp.eta_time2.b)
  }
  
  
  log.lamb <- log(lamb.old[Index0])
  log.lamb[is.na(log.lamb)] <- 0
  
  Wtime_phi.old <- if (ncw > 0) Wtime %*% phi.old else rep(0, n)
  Wtime2_phi.old <- if (ncw > 0) Wtime2 %*% phi.old else rep(0, M)
  
  
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
    
    
    eta_Btime.b_list <- lapply(1:length(Btime.b), function(i) Btime.b[[i]] * eta.old[i])
    eta_Btime.b <- Reduce('+', eta_Btime.b_list)
    
    eta_Btime2.b_list <- lapply(1:length(Btime2.b), function(i) Btime2.b[[i]] * eta.old[i])
    eta_Btime2.b <- Reduce('+', eta_Btime2.b_list)
    
    # log.density1 <- log.lamb + as.vector(Wtime_phi.old) + eta.old * Ztime.b # n*GQ matrix #
    log.density1 <- log.lamb + as.vector(Wtime_phi.old) + eta_Btime.b # n*GQ matrix #
    h.s <- as.vector(Wtime2_phi.old) + eta_Btime2.b # M*GQ matrix #  
  } else if (model == "Fixed"){
    # nB <- 0
    log.density1 <- log.lamb + as.vector(Wtime_phi.old) + eta_time.b # n*GQ matrix #
    h.s <- as.vector(Wtime2_phi.old) + eta_time2.b # M*GQ matrix #
  }else{
    stop("Invalid model type")
  } 
  
  const <- matrix(0, n, GQ) # n*GQ matrix #
  
  calc_expM2(h.s) #A = exp(A)
  temp0a <- h.s * lamb.old[Index1]
  
  const[nk != 0, ] <- calc_rowsum((Index), temp0a)
  log.density2 <- - log(1 + rho * const) # n*GQ matrix # 
  log.survival <- if (rho > 0) log.density2 / rho else - const # n*GQ matrix # 
  
  f.surv <- exp(d * log.density1 + d * log.density2 + log.survival) # n*GQ matrix #
  deno <- as.vector(f.surv %*% wGQ) # vector of length n #
  Integral <- f.surv / deno # n*GQ matrix f(bi|Oi) #
  
  f.long <- sapply(1 : n, function(i) calc_MVND(Y.st[[i]], as.vector(X.st[[i]] %*% beta.old), VY[[i]]))
  
  lgLik <- sum(log(f.long * deno / (pi ^ (ncz / 2)))) # pi^ related to sigma^2
  
  CondExp <- (1 + d * rho) / (1 + rho * const) # conditional expectation E(xi|bi,Oi), n*GQ matrix #
  
  # post.bi <- Integral %*% (t(bi) * wGQ) # n*(n*ncz) matrix #
  post.bi <- calc_M1timesM2v(Integral, bi, wGQ)
  post.bi <- if(ncz > 1) {
    t(sapply(1 : n, function(i) post.bi[i, ((i - 1) * ncz + 1) : (i * ncz)])) 
  } else { 
    matrix(diag(post.bi), nrow = n) # n*ncz matrix Ehat(bi) #
  }  
  #========== Update Bsigma ==========#
  if (ncz > 1) {
    tempB <-  fast_rbind_lapply_outerprod( bi.st )    # (n*ncz^2)*GQ matrix #      
  } else {
    tempB <- bi ^ 2
  }
  # post.bi2 <- Integral %*% (t(tempB) * wGQ) # n*(n*ncz^2) matrix #
  post.bi2 <- calc_M1timesM2v(Integral, tempB, wGQ)
  post.bi2 <- if(ncz > 1) {
    t(sapply(1 : n, function(i) post.bi2[i, ((i - 1) * ncz2 + 1) : (i * ncz2)])) 
  } else { 
    matrix(diag(post.bi2), nrow = n) # n*(ncz^2) matrix Ehat(bibi^T) #
  }
  Bsigma.new <- if (ncz > 1) matrix(colMeans(post.bi2), ncz, ncz) else mean(post.bi2) # ncz*ncz matrix #
  
  #========== Update beta: the linear regresion coefficents of regression Yi-E(Zi*bi) on X_i ==========#
  tempX <- Y - rowSums(Z * post.bi[ID, ]) # vector of length N #
  beta.old <- as.vector(qr.solve(X, tempX)); # vector of length ncx #
  beta.new <- beta.old
  
  
  #========== Update Ysigma ==========# 
  Ymu <- as.vector(X %*% beta.old) + do.call(rbind, lapply(1:n, function(i) Z.st[[i]] %*% bi.st[[i]])) # N*GQ matrix #
  
  post.resid <- ((Y - Ymu) ^ 2 * Integral[ID, ]) %*% wGQ # vector of length N #
  Ysigma2.new <- sum(post.resid) / N
  
  #========== calculate the score and gradient of phi and eta ==========#
  CondExp2 <- CondExp[nk != 0, ]
  if (model == "TimeVar"){
    temp0b <- list()
    temp1 <- list()
    temp2 <- list()
    # temp0c <- list()
    
    for(i in 1:(nB*nb)){
      temp0b[[i]] <- Btime2.b[[i]] * temp0a 
      temp1[[i]] <- CondExp2 * calc_rowsum( (Index), temp0b[[i]]) 
      for(j in 1:(nB*nb)){
        temp2[[((i-1)*(nB*nb)+j)]] <- calc_mult_rowsum2(v = Index, A = CondExp2, M = temp0b[[i]], L = Btime2.b[[j]])
      }
    }
    temp0c <- temp0b
    if (ncw > 0) {
      temp3 <- lapply(1:(ncw), function(i) calc_mult_rowsum1(v = Index, u = Wtime2[, i], A = CondExp2,  M =temp0a)) #for more than 2, they saved as list
      temp4 <- lapply(1:(ncw^2), function(i) calc_mult_rowsum1(v = Index, u = Wtime22[, i], A = CondExp2, M = temp0a)) 
    }
    
  }else if(model == "Fixed"){
    temp0b <- list()
    temp1 <- list()
    temp2 <- list()
    # temp0c <- list()
    
    for(i in 1:(nb)){
      temp0b[[i]] <- time2.b[[i]] * temp0a 
      temp1[[i]] <- CondExp2 * calc_rowsum( (Index), temp0b[[i]])
      if(nb > 1){
        for(j in 1:nb){
          temp2[[((i-1)*nb+j)]] <- calc_mult_rowsum2(v = Index, A = CondExp2, M = temp0b[[i]], L = time2.b[[j]])
        }
      }else{
        temp2[[i]] <- calc_mult_rowsum2(v = Index, A = CondExp2, M = temp0b[[i]], L = time2.b[[i]])
      }
      
    }
    temp0c <- temp0b
    if (ncw > 0) {
      temp3 <- lapply(1:(ncw), function(i) calc_mult_rowsum1(v = Index, u = Wtime2[, i], A = CondExp2,  M =temp0a)) #for more than 2, they saved as list
      temp4 <- lapply(1:(ncw^2), function(i) calc_mult_rowsum1(v = Index, u = Wtime22[, i], A = CondExp2, M = temp0a)) 
    }
    
  }
  
  if (ncw > 0){
    temp5 <- list()
    for(j in 1:length(temp0c)){
      for(k in 1:ncw){
        temp5_sub <- lapply(1:(ncw), function(i) calc_mult_rowsum1(v = Index, u = Wtime2[, i], A = CondExp2, M = temp0c[[j]]))
        temp5[[((j-1)*ncw+k)]] <- temp5_sub[[k]]
      }
    }
    
  } 
  
  Integral2 <- Integral[nk != 0,]
  post1 <- unlist(lapply(temp1, function(x) sum((x * Integral2) %*% wGQ))) # vector of length ncw #
  post2 <- unlist(lapply(temp2, function(x) sum((x * Integral2) %*% wGQ))) # vector of length ncw #
  if (ncw > 0) {
    post3 <- unlist(lapply(temp3, function(x) sum((x * Integral2) %*% wGQ))) # vector of length ncw #
    post4 <- unlist(lapply(temp4, function(x) sum((x * Integral2) %*% wGQ))) # vector of length ncw^2 #
    post5 <- unlist(lapply(temp5, function(x) sum((x * Integral2) %*% wGQ))) # vector of length ncw # 
    phiScore <- colSums(d * Wtime) - post3 # vector of length ncw #
  }
  

if (model == "TimeVar"){
  alphaScore <- numeric(nb*nB)
  for(i in 1:nB){
    for(k in 1:nb){
      ind <- nb*(i-1)+k
      alphaScore[ind] <- sum(d * (B[,i] * post.bi[,k])) - post1[ind] # Score function two terms.
    }
  }
  
}else if(model == "Fixed"){
  alphaScore <- numeric(nb)
  for(k in 1:nb){
    alphaScore[k] <- sum(d * (post.bi[,k])) - post1[k] # Score function two terms.
  }
}

if (ncw > 0) {
  nc_alpha <- length(eta.old)
  pa.score <- c(phiScore, alphaScore)
  pa.info <- matrix(0, (ncw +nc_alpha), (ncw +nc_alpha)) # (ncw+1)*(ncw+1) matrix # change dimension
  pa.info[1:ncw, 1:ncw] <- -post4 # second derivative #W part
  pa.info[((ncw + 1):(ncw+nc_alpha)), 1:ncw] <- -t(matrix(post5, ncol = nc_alpha)) # need to change cross derivative
  pa.info[1:ncw, ((ncw + 1):(ncw+nc_alpha))] <- -matrix(post5, ncol = nc_alpha)
  pa.info[((ncw + 1):(ncw+nc_alpha)), ((ncw + 1):(ncw+nc_alpha))] <- -post2 # need to change eta
  
  #=============== Update phi and alpha ===============#
  pa.old <- c(phi.old, eta.old) # vector of length (ncw+1) #
  paSVD <- svd(pa.info)
  pa.info.inv <- paSVD$v %*% diag(1/paSVD$d) %*% t(paSVD$u)
  pa.new <- pa.old - pa.info.inv %*% pa.score # vector of length (ncw+1) #
  phi.new <- pa.new[1:ncw]
  eta.new <- pa.new[(ncw + 1):(ncw+nc_alpha)]
}else {
  if(model == "TimeVar"){
    if(nb*nB > 1){
      nc_alpha <- length(eta.old)
      pa.info <- matrix(0, (nc_alpha), (nc_alpha))
      pa.info[c(1:nc_alpha), c(1:nc_alpha)] <- -post2
      #=============== Update phi and alpha ===============#
      paSVD <- svd(pa.info)
      pa.info.inv <- paSVD$v %*% diag(1/paSVD$d) %*% t(paSVD$u)
      eta.new <- eta.old - pa.info.inv %*% alphaScore 
    }else{
      eta.new <- eta.old - alphaScore / (-post2)
    }
    phi.new <- phi.old
  }else if(model == "Fixed"){
    if(nb> 1){
      nc_alpha <- length(eta.old)
      pa.info <- matrix(0, (nc_alpha), (nc_alpha))
      pa.info[c(1:nc_alpha), c(1:nc_alpha)] <- -post2
      #=============== Update phi and alpha ===============#
      paSVD <- svd(pa.info)
      pa.info.inv <- paSVD$v %*% diag(1/paSVD$d) %*% t(paSVD$u)
      eta.new <- eta.old - pa.info.inv %*% alphaScore 
    }else{
      eta.new <- eta.old - alphaScore / (-post2)
    }
    phi.new <- phi.old
  }
}

Wtime2_phi.new <- if (ncw > 0) Wtime2 %*% phi.new else rep(0, M)
#========== Calculate the new lambda with new parameters ==========# 
if (model == "TimeVar"){
  eta_Btime2.b_list_new <- lapply(1:length(Btime2.b), function(i) Btime2.b[[i]] * eta.new[i])
  eta_Btime2.b_new <- Reduce('+', eta_Btime2.b_list_new)
  h.sn <- as.vector(Wtime2_phi.new) + eta_Btime2.b_new # M*GQ matrix # 
}else if(model == "Fixed"){
  temp.eta_time2.b_new <- lapply(1:nb, function(i) time2.b[[i]]*eta.new[i])
  eta_time2.b_new <- Reduce('+', temp.eta_time2.b_new)
  h.sn <- as.vector(Wtime2_phi.new) + eta_time2.b_new # M*GQ matrix # 
}
calc_expM2(h.sn);
calc_M1_M2_M3_Hadamard(h.sn, CondExp ,  Integral, as.integer(Index - 1))
tempLamb <- calc_M_v(v = wGQ, M = h.sn)
postLamb <- calc_tapply_vect_sum(v1 = tempLamb, v2 = as.integer(Index1 - 1)); ## Check this!

lamb.new <- Index2 / postLamb

result <- list(beta = beta.new, Ysigma = sqrt(Ysigma2.new), Bsigma = Bsigma.new, phi = phi.new, 
               eta = eta.new, lamb = lamb.new, lgLik = lgLik, est.bi = post.bi, h.sn = h.sn[,1])
return(result)
} 

