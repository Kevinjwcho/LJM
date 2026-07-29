
#=============== Initial Value Calculation for Transformation Model II ===============#

InitVal_multiJEL_tv <- function (beta, model, n, n_LME, X, Z, Y, bBLUP, ID, Xtime, Ztime, Xtime2, Ztime2, Indcs, start, event, stop, W , ncw, Wtime2, rho, nk, Wtime22, d, Wtime, ID_surv, Bs) {
  
  # cvals = cntrlLst
  Index = Indcs$Index  
  Index0 = Indcs$Index0
  Index1 = Indcs$Index1
  Index2 = Indcs$Index2
  
  # tol.P = cvals$tol.P;
  # iter = cvals$max.iter;
  
  M <- lapply(Xtime2, nrow)
  nrow_surv <- nrow(W)
  
  ## Specifying longitudinal process.
  if(model == "TimeVar"){
    ## B-spline construct
    B_list = list()
    for(j in 1:n_LME){
      Bs_config <- Bs[[j]]
      Null_check = Bs_config %>% unlist
      if(is.null(Null_check))
        B_list[[j]] <- list(NULL)
      else{
        B_list[[j]] <- bs(stop, df = Bs_config$df, knots = Bs_config$knots, degree = Bs_config$degree, intercept = TRUE) %>% as.matrix
      }
    }
    # n_eta <- (ncol(B) * ncol(bBLUP))
    n_eta <- sapply(1:n_LME, function(x){
      Null_check = B_list[[x]] %>% unlist
      if(is.null(Null_check)){
        return(ncol(bBLUP[[x]]))
      }else{
        return(ncol(B_list[[x]])*ncol(bBLUP[[x]]))
      }
    }) %>% sum
    fixedOrRand_mat <- matrix(nrow = length(ID_surv), ncol = n_eta)
    fixedOrRand_mat2 <- matrix(nrow = length(Index), ncol = n_eta)
    # I have to start from this:
    eta_n <- sapply(1:n_LME, function(x){
      Null_check = B_list[[x]] %>% unlist
      if(is.null(Null_check)){
        paste0("eta_b", x, "_", c(0:(ncol(bBLUP[[x]])-1)))
      }else{
        outer(paste0("eta_B", c(1:ncol(B_list[[x]]))), paste0("_b", x, "_", c(0:(ncol(bBLUP[[x]])-1))), paste0) %>% t() %>% as.vector
      }
    }) %>% unlist
    
    var_loci = 0
    for(k in 1:n_LME){
      Null_check = B_list[[k]] %>% unlist
      if(is.null(Null_check)){
        for(j in 1:ncol(bBLUP[[k]])){
          var_loci = var_loci + 1
          fixedOrRand_mat[, var_loci] <- bBLUP[[k]][ID_surv,j]
          fixedOrRand_mat2[, var_loci] <- bBLUP[[k]][Index,j]
        }
      }else{
        for(i in 1:ncol(B_list[[k]])){
          for(j in 1:ncol(bBLUP[[k]])){
            var_loci = var_loci+1
            fixedOrRand_mat[, var_loci] <- B_list[[k]][,i]*bBLUP[[k]][ID_surv,j]
            fixedOrRand_mat2[, var_loci] <- B_list[[k]][Index,i]*bBLUP[[k]][Index,j]
          }
        }
      }
    }
    
    
    fixedOrRand.time <- fixedOrRand_mat # vector of length of lmm model #
    fixedOrRand.time2 <- fixedOrRand_mat2 # vector of length M (surv model)# 
  }else if(model == "Fixed"){
    # n_eta <- (ncol(bBLUP))
    n_eta <- sapply(bBLUP, ncol) %>% sum
    # fixedOrRand <- bBLUP[ID, ] # vector of length of lmm model #
    # fixedOrRand <- lapply(1:n_LME, function(x) bBLUP[[x]][ID[[x]],]) # vector of length of lmm model #
    # fixedOrRand.time <- bBLUP[ID_surv, ] # vector of length n_surv#
    fixedOrRand.time <- lapply(1:n_LME, function(x) bBLUP[[x]][ID_surv,]) # vector of length of lmm model #
    fixedOrRand.time <- do.call('cbind', fixedOrRand.time)
    # fixedOrRand.time2 <- bBLUP[Index, ] # vector of length M (surv model)# 
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
  } else {
    stop("Invalid model type")
  }
  
  
  
  
  #========== first fit the Cox model ==========#
  data.init <- data.frame(start = start, stop = stop, event = event, W = W, fixedOrRand = fixedOrRand.time)
  fit <- if (ncw > 0) coxph(Surv(start, stop, event) ~ W + fixedOrRand.time, data = data.init) else coxph(Surv(start, stop, event) ~ fixedOrRand.time, data = data.init)
  # fit <- if (ncw > 0) coxph(Surv(start, stop, event) ~ W, data = data.init) else coxph(Surv(start, stop, event) ~ 1, data = data.init)
  phi.old <- if (ncw > 0) fit$coefficients[1:ncw] else numeric(0)
  
  eta.old <- fit$coefficients[(ncw+1):(ncw+n_eta)] %>% as.vector
  # if(Init.method == "one"){
  #   eta.old <- c(0, 1)
  # }
  
  
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

  result <- list(phi = phi.new, eta = eta.new, lamb = lamb.new, B_list=B_list, eta_n = eta_n, ph = fit)
  return(result)
}
