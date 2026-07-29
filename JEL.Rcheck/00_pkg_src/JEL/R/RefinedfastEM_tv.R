RefinedfastEM_tv <- function(data, ph, n_LME, ID_surv, ni_surv, X.st, Z.st, Y.st, W, Ysigma, Bsigma, beta, fitLME_list, surv.init, ni, ncw,
                             B_list, n_eta, eta_n,
                             gh.nodes = 3, collect.hist = T, max.iter = 200,
                             tol=0.01, diff.type = "abs.rel", post.process = F, verbose = F){
  # Set-up ------------------------------------------------------------------

  start.time <- proc.time()[3]
  nK = n_LME
  q <- nK * 2
  diff <- 100; b.diff <- 100; iter <- 0
  uids = ID_surv; n = length(uids)
  # Data matrices ----
  # getXi getZi are functions to make a list of matrices for each subject. And matrices are stacked by nK, which is the number of longitudinal processes
  X = lapply(1:length(ni_surv), function(name) {
    sub = sapply(X.st, function(sublist) sublist[[name]], simplify = FALSE)
    return(as.matrix(Matrix::bdiag(sub)))
  })
  Z = lapply(1:length(ni_surv), function(name) {
    sub = sapply(Z.st, function(sublist) sublist[[name]], simplify = FALSE)
    return(as.matrix(Matrix::bdiag(sub)))
  })
  Y = lapply(1:length(ni_surv), function(name) {
    sub = sapply(Y.st, function(sublist) sublist[[name]], simplify = FALSE)
    return(unlist(sub) %>% matrix(ncol =1))
  })
  mi = lapply(1:length(ni_surv), function(name) {
    sapply(ni, function(sublist) sublist[[name]], simplify = TRUE)
  })


  Xk = lapply(1:length(ni_surv), function(name) {
    sub = sapply(X.st, function(sublist) sublist[[name]], simplify = FALSE)
    # sub = lapply(sub, function(x) {colnames(x) = colnames(Xtime[[1]]);return(x)})
    return(sub)
  })
  Zk = lapply(1:length(ni_surv), function(name) {
    sub = sapply(Z.st, function(sublist) sublist[[name]], simplify = FALSE)
    return(sub)
  })
  Yk = lapply(1:length(ni_surv), function(name) {
    sub = sapply(Y.st, function(sublist) sublist[[name]], simplify = TRUE)
    return(sub)
  })

  K <- lapply(1:length(ni_surv), function(x) matrix(W[x,], nrow = 1))

  # Initial conditions -----
  inits.long = list(var.e.init = Ysigma,
                    beta.init = unlist(beta),
                    D.init = as.matrix(Matrix::bdiag(Bsigma)),
                    long.fits = fitLME_list)
  # inits: initial values for gamma and eta (including time varying)
  # l0.init: initial value for baseline hazard
  # ph: coxph model fit

  inits.surv = list(inits = c(surv.init$phi, surv.init$eta),
                    l0.init = surv.init$lamb,
                    ph = surv.init$ph,
                    B = B_list)

  # MVLME for optimal initial conditions given observed data for longitudinal part
  mvlme.fit <- tryCatch({
    mvlme(data, Y, X, Z, Yk, Xk, Zk, mi, inits.long, nK)
  }, error = function(e) {
    message("Error in mvlme: ", e)
    return(NULL)
  })
  # mvlme(data, Y, X, Z, Yk, Xk, Zk, mi, inits.long, nK)
  # mvlme.fit <- mvlme(train_dataset$LMM_dat, Y, X, Z, Yk, Xk, Zk, mi, inits.long, nK)

  # Survival-related objects ----
  sv <- surv.mod(inits.surv$ph, data, l0.init = inits.surv$l0.init)
  ft <- sv$ft; nev <- sv$nev
  surv.ids <- sv$surv.ids; surv.times <- sv$surv.times
  Di <- sv$Di; Deltai.list <- as.list(Di)
  l0 <- sv$l0; l0i <- sv$l0i; l0i.list <- as.list(l0i); l0u <- sv$l0u
  Fi <- sv$Fi;

  KK <- sapply(1:n, function(x){  # For updates to \eta
    # x <- apply(K[[x]], 2, rep, nrow(Fu[[x]]))
    x <- apply(K[[x]], 2, rep, length(l0u[[x]]))
    if("numeric" %in% class(x)) x <- t(as.matrix(x))
    x
  }) # To use element wise multiplication in the update to \phi also for time-varying!

  # Extract initial conditions
  D <- mvlme.fit$D
  b <- lapply(mvlme.fit$b, c)
  beta <- c(mvlme.fit$beta)
  var.e <- mvlme.fit$var.e
  V <- mvlme.fit$V



  mvlme.time <- mvlme.fit$elapsed.time
  # gamma <- inits.surv$inits[(ncw+1):length(inits.surv$inits)]; gr <- rep(gamma)  # 2 for intercept and slope + proportional assoc
  eta <- inits.surv$inits[(ncw+1):length(inits.surv$inits)]; er <- eta
  phi <- inits.surv$inits[1:ncw]

  # Cast to parameter vector
  vD <- vech(D)
  names(vD) <- paste0('D[', apply(which(lower.tri(D, T), arr.ind = T), 1, paste0, collapse = ','),']')
  params <- c(vD, beta, var.e, eta, phi)

  # Collect data objects and iteration "0" history
  if(collect.hist) iter.hist = data.frame(iter = iter, t(params))

  # Gaussian Quadrature -----------------------------------------------------
  gh <- statmod::gauss.quad.prob(gh.nodes, "normal")
  v <- gh$n; w <- gh$w

  # Define indices for helping later
  b.inds <- split(seq(nK * 2), rep(1:nK, each = 2))
  beta.inds <- split(seq(length(beta)), rep(1:nK, each = 2))

  names(params) <- c(names(vD),
                     sapply(1:length(beta.inds), function(x) paste0('beta_', x, '_', 1:length(beta.inds[[x]])), simplify = T),
                     paste0('var.e_', 1:nK),
                     # names(eta),
                     eta_n,
                     paste0('phi_', names(phi)))

  # B_list by subject
  BB = list()
  for(k in 1:n_LME){
    Null_check = B_list[[k]] %>% unlist
    if(is.null(Null_check)){
      BB_sub = matrix(rep(1, length(ft)), ncol = 1)
    }else{
      BB_sub = predict(B_list[[k]], ft)
    }
    BB[[k]] = BB_sub
  }

  eta.inds <- sapply(BB, function(x) 2*ncol(x))
  BBc = do.call(cbind, BB)
  # Duplicate each column and place duplicates next to originals
  BBcd <- do.call(cbind, lapply(seq_len(ncol(BBc)), function(i) {
    cbind(BBc[, i], BBc[, i])
  }))

  BBi = lapply(surv.times, function(x) matrix(BBcd[x, ], ncol = length(eta_n)))

  dmats <- list(Y = Y, X = X, Z = Z, mi = mi, Yk = Yk, Xk = Xk, Zk = Zk,        # Longit.
                K = K, Deltai = Deltai.list, BB = BB, BBi = BBi)    # Survival



  message("Starting EM Algorithm")
  EM.time <- c()
  # EM ----------------------------------------------------------------------
  while(diff > tol & iter < max.iter){
    p1 <- proc.time()[3]

    ###################
    #E-step ========= #
    ###################

    etaBBs <- etaBBs_compute(BBi, eta, eta.inds)

    b.hat <- mapply(function(b, Y, X, Z, V, mi, K, Delta, l0i, l0u, etaBB){
      ucminf::ucminf(b, ll_tv, gradll_tv,
                     Y, X, Z, V, D, sum(mi), K, Delta, l0i, l0u,
                     beta, phi, etaBB, nK, q,
                     control = list(xtol = 1e-3, grtol = 1e-6))$par
    },b = b, Y = Y, X = X, Z = Z, V = V, mi = mi, K = K, Delta = Deltai.list, l0i = l0i.list,
    l0u = l0u, etaBB = etaBBs, SIMPLIFY = F)
    b.hat.split <- lapply(b.hat, function(y) lapply(b.inds, function(x) y[x]))

    Sigmai <- mapply(function(b, Z, V, K, l0u, etaBB){
      # solve(-1 * sdll_tv(b, Z, D, V, K, l0u, phi, etaBB, nK))
      MASS::ginv(-1 * sdll_tv(b, Z, D, V, K, l0u, phi, etaBB, nK))
    }, b = b.hat, Z = Z, V = V, K = K, l0u = l0u, etaBB = etaBBs, SIMPLIFY = F)

    S <- lapply(Sigmai, function(y) lapply(b.inds, function(x) y[x,x]))   # Split out into K constituent sub-matrices along block diagonal

    # Step to update D ----
    D.newi <- mapply(function(S, b){
      S + tcrossprod(b)
    }, S = Sigmai, b = b.hat, SIMPLIFY = F)

    ## Steps to update longitudinal parameters ----

    inv_omega <- calc_inv_omega(mi, var.e)
    XometX <- calc_XometX(X, inv_omega)
    Sigma.longK <- calc_Sigma_longK(Zk, S, nK)
    mu.longK <- calc_mu_longK(Yk, Xk, Zk, b.hat, beta, beta.inds, b.inds, nK)
    beta.rhs <- calc_beta_rhs(X, Y, Z, inv_omega, b.hat)
    Eee <- calc_Ee(mu.longK, Sigma.longK, nK)

    # #####
    # Survival Parameters
    # #####

    mu_surv <- calc_mu_surv_tv(X, Y, Z, inv_omega, K, D, beta, phi, etaBBs, nK)
    Sigma2_surv <- calc_Sigma2_surv_tv(X, Z, inv_omega, D, etaBBs)

    # obtain the Expectation of Wi\phi + \etab_i using weighted f(eta%*%b_i|O_i)
    Es_exp <- Esurv_exp_t(w, v, mu = mu_surv, variance = Sigma2_surv,
                          mu_new = mu_surv, variance_new = Sigma2_surv, l0i, l0u)

    # Set out Newton-Raphson items for update to (\eta, \phi)
    Sge <- Setaphi_t(c(eta, phi), X, Y, Z, inv_omega, K, D, beta, l0i, l0u, Di, nK, w, v, mu_surv, Sigma2_surv, BBi, eta.inds, eps = 0.0001)

    Hge <- Hetaphi_t(c(eta, phi), X, Y, Z, inv_omega, K, D, beta, l0i, l0u, Di, nK, w, v, mu_surv, Sigma2_surv, BBi, eta.inds, eps = 0.0001)

    # ##################
    # M-step ========= #
    # ##################

    # D -----
    D.new <- Reduce('+', D.newi)/n

    # \beta ----
    beta.new <- solve(Reduce('+', XometX)) %*% Reduce('+', beta.rhs) # NB this slightly faster than Reduce('+',.) on rhs.
    # var.e ----
    var.e.new <- colSums(do.call(rbind, Eee))/colSums(do.call(rbind, mi))      # NB this same speed as storing Ee directly as an array.
    # The baseline hazard, \lambda ----
    lambda <- update_lambda_tv(Es_exp, l0, n)

    l0.new <- nev/rowSums(do.call(cbind, lambda))
    l0u.new <- lapply(l0u, function(x){
      ll <- length(x); l0.new[1:ll]
    })
    l0i.new <- c()
    l0i.new[which(Di == 0)] <- 0
    l0i.new[which(Di == 1)] <- l0.new[match(Fi[which(Di==1), 2], ft)]

    # (\eta, \phi) ----
    eta.phi.new <- c(eta, phi) - solve(Hge, Sge)

    eta.new <- eta.phi.new[1:sum(eta.inds)]
    phi.new <- eta.phi.new[(sum(eta.inds) + 1):length(eta.phi.new)]

    EM.time[iter + 1] <- proc.time()[3] - p1 # M step finishes here ---

    # Update parameters and print ----
    params.new <- c(vech(D.new), beta.new, var.e.new, eta.new, phi.new); names(params.new) <- names(params)
    if(verbose) print(sapply(params.new, round, 4))
    # Take differences (user input)
    if(diff.type == "abs"){
      diffs <- abs(params.new-params)
      b.diff <- max(abs(do.call(rbind, b.hat) - do.call(rbind, b)))
    }else if(diff.type == "abs.rel"){
      diffs <- abs(params.new-params)/(abs(params) + 1e-3)
      b.diff <- max(abs(do.call(rbind, b.hat) - do.call(rbind, b))/(abs(do.call(rbind, b)) + 1e-3))
    }
    diff <- max(diffs)
    # Message output (max relative diff)
    message("\nIteration ", iter + 1, " maximum difference: ", round(diff, 5))
    message("Largest change: ", names(params)[which(diffs==diff)])
    message("--> old: ", params[which(diffs==diff)], " new: ", params.new[which(diffs==diff)])
    message("Largest change in random effects: ", round(b.diff, 3))

    # Update ----
    params <- params.new
    D <- D.new; var.e <- var.e.new
    V <- lapply(mi, function(iii) {
      diag(x = rep(var.e, iii), ncol = sum(iii))
    })
    eta <- eta.new
    phi <- phi.new
    beta <- beta.new;
    b <- b.hat
    l0 <- l0.new; l0u <- l0u.new; l0i <- l0i.new; l0i.list <- as.list(l0i)
    iter <- iter + 1
    if(collect.hist) iter.hist = rbind(iter.hist, c(iter = iter, t(params)))
  }
  if(collect.hist) colnames(iter.hist) = c("iter", names(params))
  # Set up and return list ----
  rownames(D) = sapply(1:length(b.inds), function(x) paste0('b', x, '_', (1:length(b.inds[[x]])-1)))
  colnames(D) = rownames(D)
  beta = as.vector(beta)
  names(beta) = sapply(1:length(beta.inds), function(x) paste0('beta', x, '_', (1:length(beta.inds[[x]])-1)), simplify = T)
  names(var.e) = paste0('var.e_', 1:nK)
  names(eta) = eta_n
  names(phi) = paste0('phi_', names(phi))
  coeffs <- list(beta = beta, var.e = var.e, D = D, eta = eta, phi = phi, hazard = cbind(ft, l0))

  rtn <- list(REs = do.call(rbind, b), coeffs = coeffs,
              dmats = dmats,
              # Elapsed times //
              EMtime = round(sum(EM.time), 2),
              mvlme.time = mvlme.fit$elapsed.time,
              comp.time = round(proc.time()[3] - start.time, 2))

  # Have to do 11/27
  if(post.process){
    message("\nStarting post-fit calculations...")
    pp.start <- proc.time()[3]

    etaBBs <- etaBBs_compute(BBi, eta, eta.inds)

    b <- mapply(function(b, Y, X, Z, V, mi, K, Delta, l0i, l0u, etaBB){
      ucminf::ucminf(b, ll_tv, gradll_tv,
                     Y, X, Z, V, D, sum(mi), K, Delta, l0i, l0u,
                     beta, phi, etaBB, nK, q,
                     control = list(xtol = 1e-3, grtol = 1e-6))$par
    },b = b, Y = Y, X = X, Z = Z, V = V, mi = mi, K = K, Delta = Deltai.list, l0i = l0i.list,
    l0u = l0u, etaBB = etaBBs, SIMPLIFY = F)

    Sigmai <- mapply(function(b, Z, V, K, l0u, etaBB){
      # solve(-1 * sdll_tv(b, Z, D, V, K, l0u, phi, etaBB, nK))
      MASS::ginv(-1 * sdll_tv(b, Z, D, V, K, l0u, phi, etaBB, nK))
      # sdll(b, Z, D, V, K, l0u, phi, etaBB, nK)
    }, b = b, Z = Z, V = V, K = K, l0u = l0u, etaBB = etaBBs, SIMPLIFY = F)
    S <- lapply(Sigmai, function(y) lapply(b.inds, function(x) y[x,x]))   # Split out into K constituent sub-matrices along block diagonal

    H <- PRES_hessian_tv(coeffs, dmats, V, b, b.hat, Sigmai, S, l0i, l0u, gh.nodes, n, q, nK, eta.inds, nev, Fi, delta = 0.0001)
    Vcov <- -solve(H)
    rownames(Vcov) = names(params)
    colnames(Vcov) = names(params)
    # SEs <- sqrt(diag(Vcov))
    # names(SEs) <- names(params)
    pp.end <- proc.time()[3]
    message("\nDone")
  }

  if(collect.hist) rtn$history <- iter.hist
  if(post.process){
    rtn$Vcov <- Vcov
    rtn$REs <- do.call(rbind, b)
    rtn$postprocess.time <- round(pp.end - pp.start, 2)
    rtn$comp.time = round(proc.time()[3] - start.time, 2)
  }
  rtn
}
