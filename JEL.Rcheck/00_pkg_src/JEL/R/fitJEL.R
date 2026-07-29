# Joint Modeling Main Function with LME (linear mixed effects)

fitJEL <- function (fitLME, fitCOX, surv_data, df = NULL, degree = 3, knots = NULL, Bknots = NULL, model = "Fixed", rho = 0, timeVarY = NULL, timeVarT = NULL,
                 Init.method = "base", control = list(), ...)
{
  # cat("Running jmodelTM(), may take some time to finish.\n")
  # call <- match.call()
  #
  # CheckInputs(fitLME, fitCOX, rho)
  #

# Sub LMM data set --------------------------------------------------------

  ID1 <- as.vector(unclass(fitLME$groups[[1]])) # Subjects
  uniqueID <- !duplicated(ID1) #
  tempID <- which(uniqueID)
  tempID <- c(tempID, length(ID1) + 1)
  ni <- diff(tempID) # the number of obs of each subjects
  ID <- rep(1:sum(uniqueID), times = ni) # ID variable
  bBLUP1 <- data.matrix(ranef(fitLME))  # random effect of each observation (slope and intercept)
  # I think this is to reconstruct the b_i random effect matrix
  if (ncol(bBLUP1) == 1) {
    bBLUP <- matrix(bBLUP1[ID1[uniqueID], ], ncol = 1)
  } else {
    bBLUP <- bBLUP1[ID1[uniqueID], ]
    dimnames(bBLUP) <- NULL
  }
  nLong <- nrow(bBLUP) # the number of observation (random effects)
  if (ncol(fitCOX$y) != 3)
    stop("\n must fit time-dependent Cox model in coxph().")



# Sub survival data set ---------------------------------------------------

  start <- as.vector(fitCOX$y[, 1])
  stop <- as.vector(fitCOX$y[, 2])
  event <- as.vector(fitCOX$y[, 3])
  ft <- stop[which(event == 1)] %>% unique() # event time

  # ID1_surv <- as.vector(unclass(fitCOX$groups[[1]])) # Subjects
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
  if (nLong != nSurv)
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
  if (!is.null(timeVarT)) {
    if (!all(timeVarT %in% all.vars(TermsSurv)))
      stop("\n'timeVarT' does not correspond columns in the fixed-effect design matrix of 'fitCOX'.")
    mfSurv[timeVarT] <- Time
  }
  if (ncw > 0) {
    Wtime <- as.matrix(model.matrix(formSurv, mfSurv)) # matrix of W with intercept
    if(attr(TermsSurv, 'intercept')) Wtime <- as.matrix(Wtime[, - 1])
    # design matrix in survival part, one row for each subject, excluding intercept #
  } else Wtime <- matrix(, ncol = 0, nrow = nSurv) # If there is no baseline variable,

  ## Sub LMM
  TermsLongX <- fitLME$terms
  mydata <- fitLME$data[all.vars(TermsLongX)] # Extract data in LMM
  formLongX <- formula(fitLME)
  mfLongX <- model.frame(TermsLongX, data = mydata)
  X <- as.matrix(model.matrix(formLongX, mfLongX)) # model matrix in LMM
  # varNames$eta.name <- rownames(attr(TermsLongX, "factors"))[attr(TermsLongX, "response")]

  formLongZ <- formula(fitLME$modelStruct$reStruct[[1]]) # random effect formula
  mfLongZ <- model.frame(terms(formLongZ), data = mydata) # Z
  TermsLongZ <- attr(mfLongZ, "terms")
  Z <- as.matrix(model.matrix(formLongZ, mfLongZ))
  Y <- as.vector(model.response(mfLongX, "numeric"))
  # give the column in mfLongX which is considered as response, may be transformed #

  data.id <- mydata[uniqueID, ] # pick the first row of each subject in mydata, nrow=n # but there is no observation time

  if (!is.null(timeVarY)) {
    if (!all(timeVarY %in% names(mydata)))
      stop("\n'timeVarY' does not correspond to columns in the fixed-effect design matrix of 'fitLME'.")
    data.id[timeVarY] <- Time
  }

  mfLongX.id <- model.frame(TermsLongX, data = data.id)
  Xtime <- as.matrix(model.matrix(formLongX, mfLongX.id)) # same structure with X, but with only n rows # unique subject

  mfLongZ.id <- model.frame(TermsLongZ, data = data.id)
  Ztime <- as.matrix(model.matrix(formLongZ, mfLongZ.id)) # same structure with Z, but with only n rows #

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
  data.id2 <- data.id[Indcs$Index, ]
  if (!is.null(timeVarY)) {
    data.id2[timeVarY] <- times
  }
  mfLongX2 <- model.frame(TermsLongX, data = data.id2)
  Xtime2 <- as.matrix(model.matrix(formLongX, mfLongX2)) # longer version than Xtime (nk)
  mfLongZ2 <- model.frame(TermsLongZ, data = data.id2)
  Ztime2 <- as.matrix(model.matrix(formLongZ, mfLongZ2)) # longer version than Ztime (nk)


  mfSurv2 <- mfSurv[Indcs$Index, ] # longer version of surv data
  if (!is.null(timeVarT)) {
    mfSurv2[timeVarT] <- times
  }
  if (ncw > 0) {
    Wtime2 <- as.matrix(model.matrix(formSurv, mfSurv2))
    if(attr(TermsSurv, 'intercept')) Wtime2 <- as.matrix(Wtime2[, - 1]) # excluding intercept #
  } else Wtime2 <- matrix(, ncol = 0, nrow = M)
  #########

  n <- nLong
  N <- length(Y)
  nu <- length(U) # number of uncensored event.
  ncz <- ncol(Z)
  ncx <- ncol(X)
  ncz2 <- ncz ^ 2
  p <- ncz * (ncz + 1) / 2 # number of Bsigma

  cntrlLst <- GenerateControlList(control, ncz)

  GHQ <- gauss.quad(cntrlLst$nknot, kind = "hermite")
  b <- as.matrix(expand.grid(rep(list(GHQ$nodes), ncz)))
  wGQ <- as.matrix(expand.grid(rep(list(GHQ$weights), ncz)))
  wGQ <- apply(wGQ, 1, prod)
  GQ <- nrow(b)

  Z.st <- lapply(split(Z, ID), function(x) matrix(x, ncol = ncz)) # Something list
  Y.st <- split(Y, ID)
  X.st <- lapply(split(X, ID), function(x) matrix(x, ncol = ncx))
  Ztime2.st <- vector('list', n)

  #####
  for (i in (1:n)[nk != 0]) { Ztime2.st[[i]] <- matrix(Ztime2[Indcs$Index == i, ], ncol = ncz) }
  #####

  Wtime22 <- if(ncw > 1) t(apply(Wtime2, 1, function(x) tcrossprod(x))) else Wtime2 ^ 2
  Xtime22 <- if(ncx > 1) t(apply(Xtime2, 1, function(x) tcrossprod(x))) else Xtime2 ^ 2
  X2 <- if(ncx > 1) t(apply(X, 1, function(x) tcrossprod(x))) else X ^ 2
  X2.sum <- matrix(colSums(X2), nrow = ncx)

  Bsigma <- lapply(lapply(fitLME$modelStruct$reStruct, as.matrix),
                   function(x) x * fitLME$sigma ^ 2)[[1]]

  # the estimated variance-covariance matrix for the random effects  #
  beta <- as.vector(fixef(fitLME))
  varNames$beta.names <- names(fixef(fitLME))
  Ysigma <- fitLME$sigma

  surv.init <- InitValJEL(beta, model = model, n = n, X = X, Z = Z, Y = Y, bBLUP = bBLUP, ID = ID, Xtime = Xtime, Ztime = Ztime, Xtime2 = Xtime2, Ztime2 = Ztime2, Indcs = Indcs,
                          start = start, event = event, stop = stop, W = W, ncw = ncw, Wtime2 = Wtime2, rho= rho, nk = nk, Wtime22 = Wtime22, d = d,  Wtime = Wtime, cvals = cntrlLst,
                          ID_surv = ID_surv, df =df, knots = knots, degree = degree, Bknots = Bknots, Init.method = Init.method )
  phi <- surv.init$phi
  eta <- surv.init$eta
  lamb <- surv.init$lamb
  B <- surv.init$B
  varNames$eta.name <- surv.init$eta_n


  theta.old <- list(beta = beta, phi = phi, eta = eta, Ysigma = Ysigma, Bsigma = Bsigma,
                    lamb = lamb, lgLik = 0)
  err.P <- err.L <- step <- 1

  err.P_list <- c(err.P)
  err.L_list <- c(err.L)
  h.sn_list <- c(0)
  while (step <= cntrlLst$max.iter) {
    # while (step <= 9) {

    if (err.P < cntrlLst$tol.P | err.L < cntrlLst$tol.L) break

    theta.new <- EMiterJEL(theta.old, n = n, Z.st = Z.st, Ztime = Ztime, Ztime2.st = Ztime2.st, nk = nk, Indcs = Indcs, Wtime2 = Wtime2, Xtime2 = Xtime2, GQ = GQ,
                           rho = rho, wGQ = wGQ, d = d, Y.st = Y.st, X.st = X.st, ncz = ncz, ncz2 = ncz2, b = b, model =  model, Wtime = Wtime, Xtime = Xtime, X = X,
                           Y = Y, ID = ID, N = N, ncw = ncw, Wtime22 = Wtime22, ncx = ncx, Xtime22 = Xtime22, Z = Z, X2.sum = X2.sum, B = B, cntrlLst = cntrlLst)
    new.P <- c(theta.new$beta, theta.new$phi, theta.new$eta, theta.new$Ysigma, theta.new$Bsigma)
    old.P <- c(theta.old$beta, theta.old$phi, theta.old$eta, theta.old$Ysigma, theta.old$Bsigma)
    h.sn_list <- c(h.sn_list, theta.new$h.sn)

    err.P <- max(abs(new.P - old.P) / (abs(old.P) + .Machine$double.eps * 2))
    err.P_list <- c(err.P_list, err.P)

    new.L <- theta.new$lgLik
    old.L <- theta.old$lgLik
    err.L <- abs(new.L - old.L) / (abs(old.L) + .Machine$double.eps * 2)
    err.L_list <- c(err.L_list, err.L)

    step <- step + 1
    theta.old <- theta.new
    cat("step:", step, "\n")
    gc()
  }
  converge <- as.numeric(err.P < cntrlLst$tol.P | err.L < cntrlLst$tol.L)

  if (cntrlLst$SE.method == 'PFDS') {
    if (CheckDeltaFD(theta.new, ncz, cntrlLst$delta)) {
      time.SE <- system.time(Vcov <- PFDS_JEL(model, theta.new, ncx = ncx, ncz = ncz, ncw = ncw, p = p, cvals = cntrlLst, varNames = varNames, Indcs = Indcs, n = n, Z.st = Z.st,
                                          Y.st = Y.st, X.st = X.st, Ztime = Ztime, nk = nk, Wtime = Wtime, Wtime2 = Wtime2, Xtime = Xtime, Xtime2 = Xtime2, GQ = GQ, rho = rho,
                                          d = d, wGQ = wGQ, ncz2 = ncz2, b = b, Ztime2.st = Ztime2.st, X = X, Y = Y, ID = ID, N = N, Z = Z, B=B))[3]
      if (any(is.na(suppressWarnings(sqrt(diag(Vcov))))))
        warning("NA's present in StdErr estimation due to numerical error!\n")
    } else {
      Vcov <- time.SE <- NA
      warning("\n 'delta' is too large, use smaller 'delta'!")
    }
  } else if (cntrlLst$SE.method == 'PRES') {
    # if (CheckDeltaRE(theta.new, ncz, cntrlLst$delta)) {
      time.SE <- system.time(Vcov <- PRES_JEL(model, theta.new, ncz = ncz, ncx = ncx, ncw = ncw, n = n, Z.st = Z.st, Y.st = Y.st, X.st = X.st, b = b, Ztime = Ztime, Ztime2.st = Ztime2.st,
                                          nk = nk, Wtime = Wtime, Xtime = Xtime, Wtime2 = Wtime2, Xtime2 = Xtime2, rho = rho, Indcs = Indcs, wGQ =wGQ, GQ = GQ, d = d, p = p, ncz2 = ncz2,
                                          X = X, Y = Y, Z = Z, ID = ID, N = N, B = B, cvals = cntrlLst, varNames = varNames))[3]
      if (any(is.na(suppressWarnings(sqrt(diag(Vcov))))))
        warning("NA's present in StdErr estimation due to numerical error!\n")
    # } else {
    #   Vcov <- time.SE <- NA
    #   warning("\n 'delta' is too large, use smaller 'delta'!")
    # }
  } else if (cntrlLst$SE.method == 'PLFD') {
    if (CheckDeltaFD(theta.new, ncz, cntrlLst$delta)) {
      time.SE <- system.time(Vcov <- PLFD_JEL(model, theta.new, n= n, ncx = ncx, ncz = ncz, ncw = ncw, p = p, cvals = cntrlLst, varNames = varNames, Z.st = Z.st, Y.st =Y.st, X.st = X.st,
                                          b = b, Ztime = Ztime, nk = nk, Indcs = Indcs, Wtime = Wtime, Xtime = Xtime, Wtime2 = Wtime2, Xtime2 = Xtime2, GQ = GQ, rho = rho, d = d,
                                          wGQ = wGQ, Ztime2.st = Ztime2.st, B = B))[3]
      if (any(is.na(suppressWarnings(sqrt(diag(Vcov))))))
        warning("NA's present in StdErr estimation due to numerical error!\n")
    } else {
      Vcov <- time.SE <- NA
      warning("\n 'delta' is too large, use smaller 'delta'!")
    }
  } else {
    Vcov <- time.SE <- NA
  warning("\n Standard error estimation method should be either 'PFDS', 'PRES' or 'PLFD'.")
  }

  theta.new$lamb <- data.frame("time" = U, "bashaz" = theta.new$lamb)
  names(theta.new$beta) <- varNames$beta.names
  names(theta.new$phi) <- varNames$phi.names
  names(theta.new$eta) <- surv.init$eta_n
  names(theta.new$Ysigma) <- "sigma.e"
  if (ncz > 1) dimnames(theta.new$Bsigma) <- dimnames(Bsigma)
  else names(theta.new$Bsigma) <- "sigma.b"
  colnames(theta.new$est.bi) <- colnames(Bsigma)
  rownames(theta.new$est.bi) <- (fitLME$groups[[1]])[uniqueID]
  if(model == "TimeVar") B <- predict(B, ft)
  else B <- NULL

  result <- list()
  result$coefficients <- theta.new
  result$logLik <- theta.new$lgLik
  result$call <- call
  result$numIter <- step
  result$Vcov <- Vcov
  result$est.bi <- theta.new$est.bi
  result$coefficients$est.bi <- NULL
  result$convergence <- if(converge == 1) "success" else "failure"
  result$control <- cntrlLst
  result$time.SE <- time.SE
  result$N <- N
  result$n <- n
  result$d <- d
  result$rho <- rho
  result$dataMat <- list(B = B, ID = ID)
  class(result) <-  unlist(strsplit(deparse(sys.call()), split = '\\('))[1]

  return(result)
}
