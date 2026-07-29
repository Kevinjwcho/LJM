# Joint Modeling Main Function with LME (linear mixed effects)

fitMultiJEL_tv <- function (fitLME_list, fitCOX, train_dataset,
                         Bs = list(
                           list(df = 3, degree = 1, knots = c(2), Bknots = NULL),
                           list(NULL),
                           list(NULL)
                         ),
                         Vcov = T, ...)
{
  # Check if fitLME_list is a list
  if (!is.list(fitLME_list)) {
    stop("fitLME_list must be a list of fitLME objects")
  }

  # Check if Bs is NULL
  if (is.null(Bs)) {
    model = "Fixed"
  }else{
    if(length(Bs) != length(fitLME_list)){
      stop("Bs must be a list of length equal to the number of LME models")
    }else{
      # if(EM.method != "RefinedfastEM"){
      #   stop("Bs is only available for RefinedfastEM method")
      # }else{
        model = "TimeVar"
      # }
    }
  }

  rho = 0 # Right now the fixed rho

  # Sub LMM data set --------------------------------------------------------
  surv_data <- train_dataset$Surv_dat
  n_LME <- length(fitLME_list)

  ID1 <- lapply(fitLME_list, function(x) as.vector(unclass(x$groups[[1]]))) # Subjects
  uniqueID <- lapply(ID1, function(x) !duplicated(x))
  tempID <- lapply(uniqueID, function(x){ y = which(x);z= c(y, length(x) + 1); return(z)})
  ni <- lapply(tempID, function(x) diff(x)) # the number of obs of each subjects
  ID <- lapply(1:n_LME, function(i) rep(1:sum(uniqueID[[i]]), times = ni[[i]])) # ID variable
  bBLUP1 <- lapply(fitLME_list, function(x) data.matrix(ranef(x)))  # random effect of each observation (slope and intercept)

  bBLUP <- lapply(1:n_LME, function(i){
    if(ncol(bBLUP1[[i]]) == 1){
      bBLUP <- matrix(bBLUP1[[i]][ID1[[i]][uniqueID[[i]]], ], ncol = 1)
    }else{
      bBLUP <- bBLUP1[[i]][ID1[[i]][uniqueID[[i]]], ]
      dimnames(bBLUP) <- NULL
    }
    return(bBLUP)
  })

  nLong <- lapply(bBLUP, function(x) nrow(x)) # the number of observation (random effects)
  if(length(unique(unlist(nLong)))>1){
    stop("The number of observations in the longitudinal data sets are not equal.")
  }else{
    nLong = unique(unlist(nLong))
  }

  if (ncol(fitCOX$y) != 3)
    stop("\n must fit time-dependent Cox model in coxph().")

  # Sub survival data set ---------------------------------------------------

  start <- as.vector(fitCOX$y[, 1])
  stop <- as.vector(fitCOX$y[, 2])
  event <- as.vector(fitCOX$y[, 3])

  ID1_surv <- as.vector(surv_data$id) # Subjects
  uniqueID_surv <- !duplicated(ID1_surv) #
  tempID_surv <- which(uniqueID_surv)
  tempID_surv <- c(tempID_surv, length(ID1_surv) + 1)
  ni_surv <- diff(tempID_surv)
  ID_surv <- rep(1:sum(uniqueID_surv), times = ni_surv) # ID variable

  Time <- stop[cumsum(ni_surv)] # last observation
  d <- event[cumsum(ni_surv)]# event index
  nSurv <- length(Time) # number of subjects
  if (sum(d) < 5)
    warning("\n more than 5 events are required.")
  # I want to write a code that nLong is a list and n is a numeric. If nLong's are not equal, stop
  if (length(unique(unlist(nLong))) > 1)
    stop("The number of observations in the longitudinal data sets are not equal.")

  if (all(unlist(nLong) != nSurv))
    stop("\n sample sizes in the longitudinal and event processes differ.")


  # Sub Survival ------------------------------------------------------------


  W <- as.matrix(fitCOX$x) # Baseline W matrix
  ncw <- ncol(W) # the number of column of W
  # phi setting
  varNames <- list()
  varNames$phi.names <- colnames(W)
  formSurv <- formula(fitCOX)
  TermsSurv <- fitCOX$terms # terms import
  mfSurv <- model.frame(TermsSurv, surv_data)[cumsum(ni_surv), ] # start, stop, event, baseline (W)
  if (ncw > 0) {
    Wtime <- as.matrix(model.matrix(formSurv, mfSurv)) # matrix of W with intercept
    if(attr(TermsSurv, 'intercept')) Wtime <- as.matrix(Wtime[, - 1])
  } else Wtime <- matrix(, ncol = 0, nrow = nSurv) # If there is no baseline variable,

  ## Sub LMM ---------------------------------------------------------------
  TermsLongX <- lapply(fitLME_list, function(x) x$terms)
  mydata <- lapply(fitLME_list, function(x) x$data[all.vars(x$terms)]) # Extract data in LMM
  formLongX <- lapply(fitLME_list, function(x) formula(x))
  mfLongX <- lapply(1:n_LME, function(i) model.frame(TermsLongX[[i]], data = mydata[[i]]))
  X <- lapply(1:n_LME, function(i) as.matrix(model.matrix(formLongX[[i]], mfLongX[[i]])))

  formLongZ <- lapply(fitLME_list, function(x) formula(x$modelStruct$reStruct[[1]])) # random effect formula
  mfLongZ <- lapply(1:n_LME, function(i) model.frame(terms(formLongZ[[i]]), data = mydata[[i]])) # Z
  TermsLongZ <- lapply(mfLongZ, function(x) attr(x, "terms"))
  Z <- lapply(1:n_LME, function(i) as.matrix(model.matrix(formLongZ[[i]], mfLongZ[[i]])))
  Y <- lapply(1:n_LME, function(i) as.vector(model.response(mfLongX[[i]], "numeric")))
  # give the column in mfLongX which is considered as response, may be transformed #

  data.id <- lapply(1:n_LME, function(i) mydata[[i]][uniqueID[[i]], ]) # pick the first row of each subject in mydata, nrow=n # but there is no observation time
  mfLongX.id <- lapply(1:n_LME, function(i) model.frame(TermsLongX[[i]], data = data.id[[i]]))
  Xtime <- lapply(1:n_LME, function(i) as.matrix(model.matrix(formLongX[[i]], mfLongX.id[[i]]))) # same structure with X, but with only n rows # unique subject

  mfLongZ.id <- lapply(1:n_LME, function(i) model.frame(TermsLongZ[[i]], data = data.id[[i]]))
  Ztime <- lapply(1:n_LME, function(i) as.matrix(model.matrix(formLongZ[[i]], mfLongZ.id[[i]]))) # same structure with Z, but with only n rows #

  ## Sub survival
  U <- sort(unique(Time[d == 1])) # ordered uncensored observed event time #
  tempU <- lapply(Time, function(t) U[t >= U])
  times <- unlist(tempU) # vector of length M # for each subject collect the past time points from all subjects.
  nk <- sapply(tempU, length)  # length of each element in times, vector of length n #
  M <- sum(nk) # sum of total uncensored time point in survival model

  Indcs <- list();

  Indcs$Index <- rep(1:nLong, nk) # repeat 1:n by nk, length M #
  Indcs$Index0 <- match(Time, U) # We can know the event time and censored time
  Indcs$Index1 <- unlist(lapply(nk[nk != 0], seq, from = 1)) # vector of length M #
  Indcs$Index2 <- colSums(d * outer(Time, U, "==")) # vector of length nu # find uncensored event time


  ## I don't know the belows are needed.
  ## M is not related to LMM.
  data.id2 <- lapply(1:n_LME, function(i) data.id[[i]][Indcs$Index, ])

  mfLongX2 <- lapply(1:n_LME, function(i) model.frame(TermsLongX[[i]], data = data.id2[[i]]))
  Xtime2 <- lapply(1:n_LME, function(i) as.matrix(model.matrix(formLongX[[i]], mfLongX2[[i]]))) # longer version than Xtime (nk)
  mfLongZ2 <- lapply(1:n_LME, function(i) model.frame(TermsLongZ[[i]], data = data.id2[[i]]))
  Ztime2 <- lapply(1:n_LME, function(i) as.matrix(model.matrix(formLongZ[[i]], mfLongZ2[[i]]))) # longer version than Ztime (nk)


  mfSurv2 <- mfSurv[Indcs$Index, ] # longer version of surv data
  if (ncw > 0) {
    Wtime2 <- as.matrix(model.matrix(formSurv, mfSurv2))
    if(attr(TermsSurv, 'intercept')) Wtime2 <- as.matrix(Wtime2[, - 1]) # excluding intercept #
  } else Wtime2 <- matrix(, ncol = 0, nrow = M)
  #########

  n <- nLong
  N <- sum(sapply(Y, function(x) length(x)))
  nu <- length(U) # number of uncensored event.
  ncz <- sum(sapply(Z, function(x) ncol(x)))
  ncx <- sum(sapply(X, ncol))
  ncz2 <- ncz ^ 2
  p <- ncz * (ncz + 1) / 2 # number of Bsigma

  # GHQ <- gauss.quad(cntrlLst$nknot, kind = "hermite")
  # b <- as.matrix(expand.grid(rep(list(GHQ$nodes), ncz)))
  # wGQ <- as.matrix(expand.grid(rep(list(GHQ$weights), ncz)))
  # wGQ <- apply(wGQ, 1, prod)
  # GQ <- nrow(b)

  Z.st <- lapply(1:n_LME, function(i){
    lapply(split(Z[[i]], ID[[i]]), function(x) matrix(x, ncol = ncol(Z[[i]])))
  })
  Y.st <- lapply(1:n_LME, function(i) split(Y[[i]], ID[[i]]))
  X.st <- lapply(1:n_LME, function(i){
    lapply(split(X[[i]], ID[[i]]), function(x) matrix(x, ncol = ncol(X[[i]])))
  })


  Ztime2.st <- lapply(1:n_LME, function(i){
    Ztime2.sub <- vector('list', n)
    for (ii in (1:n)[nk != 0]) { Ztime2.sub[[ii]] <- matrix(Ztime2[[i]][Indcs$Index[[i]] == ii, ], ncol = ncol(Z[[i]])) }
    return(Ztime2.sub)
  })

  Wtime22 <- if(ncw > 1) t(apply(Wtime2, 1, function(x) tcrossprod(x))) else Wtime2 ^ 2
  Xtime22 <- lapply(1:n_LME, function(i){
    if(ncx > 1) t(apply(Xtime2[[i]], 1, function(x) tcrossprod(x))) else Xtime2[[i]] ^ 2
  })
  X2 <- lapply(1:n_LME, function(i){
    if(ncx > 1) t(apply(X[[i]], 1, function(x) tcrossprod(x))) else X[[i]] ^ 2
  })
  X2.sum <- lapply(1:n_LME, function(i) matrix(colSums(X2[[i]]), nrow = ncol(X[[i]])))

  Bsigma <- lapply(fitLME_list, function(x) lapply(lapply(x$modelStruct$reStruct, as.matrix),
                                                   function(y) y * x$sigma ^ 2)[[1]])


  # the estimated variance-covariance matrix for the random effects  #
  beta <- lapply(fitLME_list, function(x) as.vector(fixef(x)))
  varNames$beta.names <- lapply(fitLME_list, function(x) names(fixef(x)))
  Ysigma <- sapply(fitLME_list, function(x) x$sigma)

  surv.init <- InitVal_multiJEL_tv(beta, model = model, n = n, n_LME = n_LME, X = X, Z = Z, Y = Y, bBLUP = bBLUP, ID = ID, Xtime = Xtime, Ztime = Ztime, Xtime2 = Xtime2, Ztime2 = Ztime2, Indcs = Indcs,
                                   start = start, event = event, stop = stop, W = W, ncw = ncw, Wtime2 = Wtime2, rho= rho, nk = nk, Wtime22 = Wtime22, d = d,  Wtime = Wtime,
                                   ID_surv = ID_surv, Bs = Bs)
  phi <- surv.init$phi
  eta <- surv.init$eta
  lamb <- surv.init$lamb
  B_list <- surv.init$B_list
  eta_n <- surv.init$eta_n
  varNames$eta.name <- surv.init$eta_n


  theta.old <- list(beta = beta, phi = phi, eta = eta, Ysigma = Ysigma, Bsigma = Bsigma,
                    lamb = lamb, lgLik = 0)
  err.P <- err.L <- step <- 1

  err.P_list <- c(err.P)
  err.L_list <- c(err.L)
  h.sn_list <- c(0)
  # if(EM.method == "RefinedfastEM"){
    theta.refinedfastEM <- RefinedfastEM_tv(surv_data, surv.init$ph, n_LME, ID_surv, ni_surv, X.st, Z.st, Y.st, W, Ysigma, Bsigma, beta, fitLME_list, surv.init, ni, ncw,
                                            B_list, n_eta, eta_n, post.process = Vcov)
    gc()

    # return(theta.fastEM)
    if(nrow(theta.refinedfastEM$history) <200){
      converge = 1
    }else{
      converge = 0
    }

    result <- list()
    lamb.new = as.data.frame(theta.refinedfastEM$coeffs$hazard)
    colnames(lamb.new) = c("time", "bashaz")
    result$coefficients <- list(beta = theta.refinedfastEM$coeffs$beta,
                                phi = theta.refinedfastEM$coeffs$phi,
                                eta = theta.refinedfastEM$coeffs$eta,
                                Ysigma = sqrt(theta.refinedfastEM$coeffs$var.e),
                                Bsigma = theta.refinedfastEM$coeffs$D,
                                lamb = lamb.new)
    result$logLik <- NULL
    result$call <- call
    # result$numIter <- step
    result$Vcov <- theta.refinedfastEM$Vcov
    result$est.bi <- theta.refinedfastEM$REs
    result$convergence <- if(converge == 1) "success" else "failure"
    # result$control <- cntrlLst
    result$time.SE <- theta.refinedfastEM$postprocess.time
    result$N <- N
    result$n <- n
    result$d <- d
    # result$rho <- rho
    result$dataMat <- list(B = theta.refinedfastEM$dmats$BB, ID = ID)
    # ,   # IDName = fitLME$groups[[1]])
    class(result) <-  unlist(strsplit(deparse(sys.call()), split = '\\('))[1]
    return(result)
}
