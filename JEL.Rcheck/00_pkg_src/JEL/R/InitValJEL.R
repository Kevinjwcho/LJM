
#=============== Initial Value Calculation for Transformation Model II ===============#

InitValJEL <- function (beta, model, n, X, Z, Y, bBLUP, ID, Xtime, Ztime, Xtime2, Ztime2, Indcs, start, event, stop, W , ncw, Wtime2, rho, nk, Wtime22, d, Wtime, cvals, ID_surv, knots, df, degree, Bknots,
                        Init.method = "base") {

  Index = Indcs$Index
  Index0 = Indcs$Index0
  Index1 = Indcs$Index1
  Index2 = Indcs$Index2

  tol.P = cvals$tol.P;
  iter = cvals$max.iter;

  M <- nrow(Xtime2)
  nrow_surv <- nrow(W)
  # b_degree = 3

  if(Init.method == "mean"){
    ## For centered observation
    ID_split <- split(1:length(ID), ID)
    last_ind <- sapply(ID_split, function(x) x[length(x)])
    meanY <- mean(Y[last_ind])
    bBLUP[,1] <- Y[last_ind] - meanY
  }


  ## Specifying longitudinal process.
  if(model == "TimeVar"){
    ## B-spline construct
    if(is.null(Bknots)){
      B <- bs(stop, df = df, knots = knots, degree = degree, intercept = TRUE) %>% as.matrix
    }else{
      B <- bs(stop, df = df, knots = knots, degree = degree, Boundary.knots = Bknots, intercept = TRUE) %>% as.matrix
    }

    # B.time2 <- bs(stop[Index]) %>% as.matrix
    n_eta <- (ncol(B) * ncol(bBLUP))
    fixedOrRand_mat <- matrix(nrow = length(ID_surv), ncol = n_eta)
    fixedOrRand_mat2 <- matrix(nrow = length(Index), ncol = n_eta)
    if(ncol(bBLUP) == 1){
      eta_n <- paste0("eta_B_", c(1:ncol(B)), "_b0")
    }else{
      eta_n <- outer(paste0("eta_B", c(1:ncol(B))), paste0("_b_", c(0:(ncol(bBLUP)-1))), paste0) %>% t() %>% as.vector
    }
    for(i in 1:ncol(B)){
      for(j in 1:ncol(bBLUP)){
        fixedOrRand_mat[, ((i-1)*ncol(bBLUP)+j)] <- B[,i]*bBLUP[ID_surv,j]
        fixedOrRand_mat2[, ((i-1)*ncol(bBLUP)+j)] <- B[Index,i]*bBLUP[Index,j]
      }
    }

    fixedOrRand.time <- fixedOrRand_mat # vector of length of lmm model #
    fixedOrRand.time2 <- fixedOrRand_mat2 # vector of length M (surv model)#
  }else if(model == "Fixed"){
      n_eta <- (ncol(bBLUP))
      fixedOrRand <- bBLUP[ID, ] # vector of length of lmm model #
      fixedOrRand.time <- bBLUP[ID_surv, ] # vector of length n_surv#
      fixedOrRand.time2 <- bBLUP[Index, ] # vector of length M (surv model)#

    if(ncol(bBLUP) == 1){
      eta_n <- paste0("eta_b0")
    }else{
      eta_n <- paste0("eta_b_", c(0:(ncol(bBLUP)-1)))  %>% as.vector
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
  if(Init.method == "one"){
    eta.old <- c(0, 1)
  }


  Wtime2_phi.old <- if (ncw > 0) Wtime2 %*% phi.old else rep(0, M)
  if(n_eta == 1){
    temp <- as.vector(exp(Wtime2_phi.old + fixedOrRand.time2*eta.old)) # M*1 vector #
  } else{
    temp <- as.vector(exp(Wtime2_phi.old + fixedOrRand.time2%*%eta.old)) # M*1 vector #
  }

  lamb.old <- Index2 / calc_tapply_vect_sum(v1 = temp, v2 = as.integer(Index1 - 1)) # vector of length n_u # Index 1 with length M.


  if (rho == 0) {
    phi.new <- phi.old
    eta.new <- eta.old
    lamb.new <- lamb.old
  } else { # we don't have to see the below, since we only consider rho = 0.
    for (it in 1:iter) {
      if(n_eta == 1){
        exp.es <- exp(as.vector(Wtime2_phi.old + fixedOrRand.time2*eta.old)) # M*1 vector #
      } else{
        exp.es <- exp(as.vector(Wtime2_phi.old + fixedOrRand.time2%*%eta.old)) # M*1 vector #
      }
      const <- rep(0, n)
      temp0a <- exp.es * lamb.old[Index1];
      const[nk != 0] <- calc_tapply_vect_sum(v1 = temp0a, v2 = as.integer(Index - 1))[nk != 0] # vector of length n #

      CondExp <- (1 + d * rho) / (1 + rho * const) # conditional expectation E(xi|Oi), vector of length n #
      CondExp2 <- CondExp[nk != 0]

      temp0b <- fixedOrRand.time2 * temp0a;



      if(n_eta == 1){
        temp1 <- sum(CondExp2 *calc_tapply_vect_sum(  v1 = temp0b, v2 = as.integer(Index - 1)))
        temp2 <- sum(CondExp2 *calc_tapply_vect_sum(  v1 = fixedOrRand.time2 * temp0b, v2 = as.integer(Index - 1)));
      }else{
        temp1 <- list()
        temp2 <- list()
        for(i in 1:n_eta){
          temp1[[i]] <- sum(CondExp2 *calc_tapply_vect_sum(  v1 = temp0b[,i], v2 = as.integer(Index - 1))[nk != 0])
          for(j in 1:n_eta){
            temp2[[((i-1)*n_eta + j)]] <- CondExp2 * calc_tapply_vect_sum(v1 = fixedOrRand.time2[,j] * temp0b[,i], v2 = as.integer(Index - 1))[nk != 0]
          }
        }
        temp2 <- sapply(temp2, sum)
      }

      if (ncw > 0) {
        temp3 <- lapply(1:ncw, function(i) CondExp2 * calc_tapply_vect_sum(  v1 = Wtime2[, i] * temp0a, v2 =  as.integer(Index - 1))[nk != 0])
        temp3 <- sapply(temp3, sum) # vector of length ncw #
        temp4 <-  lapply(1:ncw^2, function(i) CondExp2 * calc_tapply_vect_sum(  v1 = Wtime22[, i] * temp0a, v2 =  as.integer(Index - 1))[nk != 0])
        temp4 <- sapply(temp4, sum) # vector of length ncw^2 #
        if(n_eta == 1){
          temp5 <- lapply(1:ncw, function(i) CondExp2 * calc_tapply_vect_sum(  v1 = Wtime2[, i] * temp0b, v2 =  as.integer(Index - 1))[nk != 0])
          temp5 <- sapply(temp5, sum) # vector of length ncw #
        }else{
          temp5 <- list()
          for(i in 1:n_eta){
            for(j in 1:ncw){
              temp5[[((i-1)*ncw + j)]] <- CondExp2 * calc_tapply_vect_sum(  v1 = Wtime2[, j] * temp0b[,i], v2 =  as.integer(Index - 1))[nk != 0]
            }
          }
          temp5 <- sapply(temp5, sum) # vector of length ncw #
        }


        phiScore <- colSums(d * Wtime) - temp3 # vector of length ncw #
      }
      if(n_eta == 1){
        alphaScore <- sum(d * fixedOrRand.time) - temp1
      }else{
        alphaScore <- numeric(n_eta)
        for(i in 1:n_eta){
          alphaScore[i] <- sum(d * fixedOrRand.time[,i]) - temp1[[i]]
        }

      }


      if (ncw > 0) {
        pa.score <- c(phiScore, alphaScore)
        pa.info <- matrix(0, (ncw + n_eta), (ncw + n_eta)) # (ncw+n_eta)*(ncw+n_eta) matrix #
        pa.info[1:ncw, 1:ncw] <- - temp4
        if(n_eta == 1){
          pa.info[((ncw + 1):(ncw+n_eta)), 1:ncw] <- - temp5
          pa.info[1:ncw, ((ncw + 1):(ncw + n_eta))] <- - temp5
        }else{
          pa.info[((ncw + 1):(ncw+n_eta)), 1:ncw] <- - t(matrix(temp5, ncol = n_eta))
          pa.info[1:ncw, ((ncw + 1):(ncw + n_eta))] <- - (matrix(temp5, ncol = n_eta))
        }

        pa.info[((ncw + 1):(ncw + n_eta)), ((ncw + 1):(ncw + n_eta))] <- - temp2

        #=============== Update phi and alpha ===============#
        pa.old <- c(phi.old, eta.old) # vector of length (ncw+1) #
        paSVD <- svd(pa.info)
        pa.info.inv <- paSVD$v %*% diag(1/paSVD$d) %*% t(paSVD$u)
        pa.new <- pa.old - pa.info.inv %*% pa.score # vector of length (ncw+1) #
        phi.new <- pa.new[1:ncw]
        eta.new <- pa.new[(ncw + 1):(ncw + n_eta)]
      } else {
        if(n_eta == 1){
          eta.new <- eta.old - alphaScore / (-temp2)
          phi.new <- phi.old
          pa.new <- eta.new
          pa.old <- eta.old
        }else{
          phi.new <- phi.old
          pa.score <- alphaScore
          pa.old <- c(eta.old)
          pa.info <- matrix(0, (n_eta), (n_eta)) # (ncw+n_eta)*(ncw+n_eta) matrix #
          pa.info[(1:n_eta), (1:n_eta)] <- - temp2

          paSVD <- svd(pa.info)
          pa.info.inv <- paSVD$v %*% diag(1/paSVD$d) %*% t(paSVD$u)
          pa.new <- pa.old - pa.info.inv %*% pa.score # vector of length (ncw+1) #
          eta.new <- pa.new[1:n_eta]
        }
      }

      Wtime2_phi.new <- if (ncw > 0) Wtime2 %*% phi.new else rep(0, M)
      #========== Calculate the new lambda with new parameters ==========#
      if(n_eta == 1){
        exp.esn <- exp(as.vector(Wtime2_phi.new + eta.new * fixedOrRand.time2))
      } else{
        exp.esn <- exp(as.vector(Wtime2_phi.new + fixedOrRand.time2%*%eta.new))
      }

      tempLamb <- calc_tapply_vect_sum(  v1 = CondExp[Index] * exp.esn, v2 = as.integer(Index1 - 1))
      lamb.new <- Index2 / tempLamb

      #========== Check Convergence ==========#
      err <- max(abs(pa.new - pa.old) / (abs(pa.old) + tol.P))
      if (err <= tol.P) break
      else {
        phi.old <- phi.new
        eta.old <- eta.new
        lamb.old <- lamb.new
        Wtime2_phi.old <- if (ncw > 0) Wtime2 %*% phi.old else rep(0, M)
      }
    }
  }

  if(model == "Fixed") B <- NULL
  result <- list(phi = phi.new, eta = eta.new, lamb = lamb.new, B=B, eta_n = eta_n)
  return(result)
}
