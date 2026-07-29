fitLLAJEL <- function(prep, fitCOX, train_dataset,
                      model = "Fixed", Vcov = TRUE,
                      gh.nodes = 5, collect.hist = TRUE,
                      max.iter = 200, tol = 0.01, diff.type = "abs.rel",
                      post.process = Vcov, verbose = FALSE, ...) {
  
  # -----------------------------
  # Checks
  # -----------------------------
  if (!is.list(prep)) stop("prep must be a list returned by prep_LLA_landmark().")
  if (is.null(prep$surv_s2)) stop("prep$surv_s2 is missing.")
  if (is.null(prep$Y) || is.null(prep$Z) || is.null(prep$W_ker)) stop("prep must contain Y/Z/W_ker.")
  if (is.null(prep$ID) || is.null(prep$ni) || is.null(prep$uniqueID)) stop("prep must contain ID/ni/uniqueID.")
  if (is.null(prep$bLLA) || is.null(prep$cLLA)) stop("prep must contain bLLA and cLLA.")
  if (is.null(prep$subject_ids)) stop("prep$subject_ids is missing.")
  if (!inherits(fitCOX, "coxph")) stop("fitCOX must be a coxph object.")
  if (ncol(fitCOX$y) != 3) stop("fitCOX must be Surv(start, stop, event) (time-dependent Cox).")
  if (is.null(fitCOX$x)) stop("fitCOX$x is missing. Fit Cox with x=TRUE, model=TRUE.")
  if (!is.list(train_dataset) || is.null(train_dataset$LMM_dat)) stop("train_dataset$LMM_dat is required.")
  if (!exists("InitVal_LLAJEL", mode = "function")) stop("InitVal_LLAJEL() not found.")
  if (!exists("RefinedfastEM_LLA", mode = "function")) stop("RefinedfastEM_LLA() not found.")
  
  surv_data <- prep$surv_s2
  
  # -----------------------------
  # Survival: extract start/stop/event and subject indexing
  # -----------------------------
  start <- as.vector(fitCOX$y[, 1])
  stop  <- as.vector(fitCOX$y[, 2])
  event <- as.vector(fitCOX$y[, 3])
  
  ID1_surv <- as.vector(surv_data$id)
  uniqueID_surv <- !duplicated(ID1_surv)
  idx_surv <- which(uniqueID_surv)
  ni_surv <- diff(c(idx_surv, length(ID1_surv) + 1))
  ID_surv <- rep(seq_len(sum(uniqueID_surv)), times = ni_surv)
  
  Time <- stop[cumsum(ni_surv)]
  d    <- event[cumsum(ni_surv)]
  nSurv <- length(Time)
  
  if (sum(d) < 5) warning("More than 5 events are recommended.")
  
  # -----------------------------
  # Cox design: Wtime (subject-level, no intercept)  [matches fitMultiJEL]
  # -----------------------------
  formSurv  <- formula(fitCOX)
  TermsSurv <- fitCOX$terms
  mfSurv <- model.frame(TermsSurv, surv_data)[cumsum(ni_surv), , drop = FALSE]
  
  W0 <- as.matrix(fitCOX$x)
  ncw0 <- ncol(W0)
  
  if (ncw0 > 0) {
    Wtime <- as.matrix(model.matrix(formSurv, mfSurv))
    if (attr(TermsSurv, "intercept")) Wtime <- Wtime[, -1, drop = FALSE]
  } else {
    Wtime <- matrix(nrow = nSurv, ncol = 0)
  }
  ncw <- ncol(Wtime)
  
  # -----------------------------
  # Build Indcs / Wtime2 for the fixed (non-time-varying) model (no longitudinal time blocks)
  # -----------------------------
  # U: ordered uncensored event times at subject-level
  U <- sort(unique(Time[d == 1]))
  tempU <- lapply(Time, function(t) U[t >= U])
  times <- unlist(tempU)
  nk <- sapply(tempU, length)
  M <- sum(nk)
  
  Indcs <- list()
  Indcs$Index  <- rep(seq_len(nSurv), nk)             # length M
  Indcs$Index0 <- match(Time, U)                      # length nSurv
  Indcs$Index1 <- unlist(lapply(nk[nk != 0], seq, from = 1))
  Indcs$Index2 <- colSums(d * outer(Time, U, "=="))
  
  # Expanded Wtime2 (risk-set expanded)  [matches fitMultiJEL]
  mfSurv2 <- mfSurv[Indcs$Index, , drop = FALSE]
  if (ncw > 0) {
    Wtime2 <- as.matrix(model.matrix(formSurv, mfSurv2))
    if (attr(TermsSurv, "intercept")) Wtime2 <- Wtime2[, -1, drop = FALSE]
  } else {
    Wtime2 <- matrix(nrow = M, ncol = 0)
  }
  
  # -----------------------------
  # Longitudinal objects from prep -> format expected by RefinedfastEM_LLA
  # -----------------------------
  y_vars     <- prep$y_vars
  n_LLA      <- length(y_vars)
  n_subject  <- length(prep$subject_ids)
  
  if (n_subject != nSurv) {
    stop("Mismatch: n_subject from prep != nSurv from survival data (unique id count).")
  }
  
  # Bsigma: allow either block matrix (2*n_LLA x 2*n_LLA) or list of 2x2 blocks
  Bsigma_in <- prep$Bsigma
  if (is.matrix(Bsigma_in)) {
    if (nrow(Bsigma_in) != 2*n_LLA || ncol(Bsigma_in) != 2*n_LLA) {
      stop("prep$Bsigma matrix must be (2*n_LLA) x (2*n_LLA).")
    }
    Bsigma <- lapply(seq_len(n_LLA), function(k) {
      idx <- (2*k-1):(2*k)
      Bsigma_in[idx, idx, drop = FALSE]
    })
  } else if (is.list(Bsigma_in)) {
    Bsigma <- Bsigma_in
  } else {
    stop("prep$Bsigma must be a matrix or a list of 2x2 blocks.")
  }
  
  # Z.st / Y.st / Wker.st (list length n_LLA; each is list length n_subject)
  Z.st <- lapply(seq_len(n_LLA), function(k) split(prep$Z[[k]], prep$ID[[k]]))
  Z.st <- lapply(seq_len(n_LLA), function(k) {
    lapply(Z.st[[k]], function(x) matrix(x, ncol = ncol(prep$Z[[k]])))
  })
  Y.st <- lapply(seq_len(n_LLA), function(k) split(prep$Y[[k]], prep$ID[[k]]))
  Wker.st <- lapply(seq_len(n_LLA), function(k) split(prep$W_ker[[k]], prep$ID[[k]]))
  
  # Initial longitudinal pieces from prep
  Ysigma2 <- as.numeric(prep$Ysigma2)
  bLLA    <- prep$bLLA
  cLLA    <- prep$cLLA
  ni      <- prep$ni
  
  # -----------------------------
  # Build surv.init INSIDE fitLLAJEL (your request)
  # -----------------------------
  surv.init <- InitVal_LLAJEL(
    bBLUP      = bLLA,
    n_subject  = n_subject,
    n_LME      = n_LLA,
    Indcs      = Indcs,
    start      = start,
    stop       = stop,
    event      = event,
    W          = Wtime,
    ncw        = ncw,
    Wtime2     = Wtime2,
    ID_surv    = ID_surv
  )
  
  # -----------------------------
  # Run EM (LLA)
  # -----------------------------
  theta.refinedfastEM <- RefinedfastEM_LLA(
    data         = train_dataset$LMM_dat,
    n_LLA        = n_LLA,
    ID_surv      = ID_surv,
    ni_surv      = ni_surv,
    Z.st         = Z.st,
    Y.st         = Y.st,
    Wker.st      = Wker.st,
    Wtime        = Wtime,
    Ysigma2      = Ysigma2,
    Bsigma       = Bsigma,
    bLLA         = bLLA,
    cLLA         = cLLA,
    n_subject    = n_subject,
    surv.init    = surv.init,
    ni           = ni,
    ncw          = ncw,
    gh.nodes     = gh.nodes,
    collect.hist = collect.hist,
    max.iter     = max.iter,
    tol          = tol,
    diff.type    = diff.type,
    post.process = post.process,
    verbose      = verbose
  )
  gc()
  
  converge <- if (!is.null(theta.refinedfastEM$history) && nrow(theta.refinedfastEM$history) < max.iter) 1 else 0
  
  # -----------------------------
  # Output (keep JEL-style structure)
  # -----------------------------
  result <- list()
  lamb.new <- as.data.frame(theta.refinedfastEM$coeffs$hazard)
  colnames(lamb.new) <- c("time", "bashaz")
  
  result$coefficients <- list(
    phi    = theta.refinedfastEM$coeffs$phi,
    eta    = theta.refinedfastEM$coeffs$eta,
    Ysigma = sqrt(theta.refinedfastEM$coeffs$var.e),
    Bsigma = theta.refinedfastEM$coeffs$D,
    lamb   = lamb.new
  )
  
  result$logLik <- NULL
  result$call <- match.call()
  if (!is.null(theta.refinedfastEM$Vcov)) result$Vcov <- theta.refinedfastEM$Vcov
  result$est.bi <- theta.refinedfastEM$REs
  result$convergence <- if (converge == 1) "success" else "failure"
  if (!is.null(theta.refinedfastEM$postprocess.time)) result$time.SE <- theta.refinedfastEM$postprocess.time
  
  result$N <- sum(sapply(prep$Y, length))
  result$n <- n_subject
  result$d <- d
  result$dataMat <- list(ID_surv = ID_surv, Indcs = Indcs)
  
  class(result) <- "JEL"
  result
}

# Time-varying version of fitLLAJEL
fitLLAJEL_tv <- function(prep, fitCOX, train_dataset,
                         Bs = NULL,
                         Vcov = TRUE,
                         gh.nodes = 5, collect.hist = TRUE,
                         max.iter = 200, tol = 0.01, diff.type = "abs.rel",
                         post.process = Vcov, verbose = FALSE, ...) {

  # -----------------------------
  # Checks
  # -----------------------------
  if (!is.list(prep)) stop("prep must be a list returned by prep_LLA_landmark().")
  if (is.null(prep$surv_s2)) stop("prep$surv_s2 is missing.")
  if (is.null(prep$Y) || is.null(prep$Z) || is.null(prep$W_ker)) stop("prep must contain Y/Z/W_ker.")
  if (is.null(prep$ID) || is.null(prep$ni) || is.null(prep$uniqueID)) stop("prep must contain ID/ni/uniqueID.")
  if (is.null(prep$bLLA) || is.null(prep$cLLA)) stop("prep must contain bLLA and cLLA.")
  if (is.null(prep$subject_ids)) stop("prep$subject_ids is missing.")
  if (!inherits(fitCOX, "coxph")) stop("fitCOX must be a coxph object.")
  if (ncol(fitCOX$y) != 3) stop("fitCOX must be Surv(start, stop, event) (time-dependent Cox).")
  if (is.null(fitCOX$x)) stop("fitCOX$x is missing. Fit Cox with x=TRUE, model=TRUE.")
  if (!is.list(train_dataset) || is.null(train_dataset$LMM_dat)) stop("train_dataset$LMM_dat is required.")
  if (is.null(Bs)) stop("Bs must be provided for time-varying model.")

  surv_data <- prep$surv_s2

  # -----------------------------
  # Survival: extract start/stop/event and subject indexing
  # (same as fitLLAJEL)
  # -----------------------------
  start <- as.vector(fitCOX$y[, 1])
  stop  <- as.vector(fitCOX$y[, 2])
  event <- as.vector(fitCOX$y[, 3])

  ID1_surv <- as.vector(surv_data$id)
  uniqueID_surv <- !duplicated(ID1_surv)
  idx_surv <- which(uniqueID_surv)
  ni_surv <- diff(c(idx_surv, length(ID1_surv) + 1))
  ID_surv <- rep(seq_len(sum(uniqueID_surv)), times = ni_surv)

  Time <- stop[cumsum(ni_surv)]
  d    <- event[cumsum(ni_surv)]
  nSurv <- length(Time)

  if (sum(d) < 5) warning("More than 5 events are recommended.")

  # -----------------------------
  # Cox design: Wtime
  # -----------------------------
  formSurv  <- formula(fitCOX)
  TermsSurv <- fitCOX$terms
  mfSurv <- model.frame(TermsSurv, surv_data)[cumsum(ni_surv), , drop = FALSE]

  W0 <- as.matrix(fitCOX$x)
  ncw0 <- ncol(W0)

  if (ncw0 > 0) {
    Wtime <- as.matrix(model.matrix(formSurv, mfSurv))
    if (attr(TermsSurv, "intercept")) Wtime <- Wtime[, -1, drop = FALSE]
  } else {
    Wtime <- matrix(nrow = nSurv, ncol = 0)
  }
  ncw <- ncol(Wtime)

  # -----------------------------
  # Build Indcs / Wtime2
  # -----------------------------
  U <- sort(unique(Time[d == 1]))
  tempU <- lapply(Time, function(t) U[t >= U])
  times <- unlist(tempU)
  nk <- sapply(tempU, length)
  M <- sum(nk)

  Indcs <- list()
  Indcs$Index  <- rep(seq_len(nSurv), nk)
  Indcs$Index0 <- match(Time, U)
  Indcs$Index1 <- unlist(lapply(nk[nk != 0], seq, from = 1))
  Indcs$Index2 <- colSums(d * outer(Time, U, "=="))

  mfSurv2 <- mfSurv[Indcs$Index, , drop = FALSE]
  if (ncw > 0) {
    Wtime2 <- as.matrix(model.matrix(formSurv, mfSurv2))
    if (attr(TermsSurv, "intercept")) Wtime2 <- Wtime2[, -1, drop = FALSE]
  } else {
    Wtime2 <- matrix(nrow = M, ncol = 0)
  }

  # -----------------------------
  # Longitudinal objects from prep
  # -----------------------------
  y_vars     <- prep$y_vars
  n_LLA      <- length(y_vars)
  n_subject  <- length(prep$subject_ids)

  if (n_subject != nSurv) {
    stop("Mismatch: n_subject from prep != nSurv from survival data.")
  }

  Bsigma_in <- prep$Bsigma
  if (is.matrix(Bsigma_in)) {
    if (nrow(Bsigma_in) != 2*n_LLA || ncol(Bsigma_in) != 2*n_LLA) {
      stop("prep$Bsigma matrix must be (2*n_LLA) x (2*n_LLA).")
    }
    Bsigma <- lapply(seq_len(n_LLA), function(k) {
      idx <- (2*k-1):(2*k)
      Bsigma_in[idx, idx, drop = FALSE]
    })
  } else if (is.list(Bsigma_in)) {
    Bsigma <- Bsigma_in
  } else {
    stop("prep$Bsigma must be a matrix or a list of 2x2 blocks.")
  }

  Z.st <- lapply(seq_len(n_LLA), function(k) split(prep$Z[[k]], prep$ID[[k]]))
  Z.st <- lapply(seq_len(n_LLA), function(k) {
    lapply(Z.st[[k]], function(x) matrix(x, ncol = ncol(prep$Z[[k]])))
  })
  Y.st <- lapply(seq_len(n_LLA), function(k) split(prep$Y[[k]], prep$ID[[k]]))
  Wker.st <- lapply(seq_len(n_LLA), function(k) split(prep$W_ker[[k]], prep$ID[[k]]))

  Ysigma2 <- as.numeric(prep$Ysigma2)
  bLLA    <- prep$bLLA
  cLLA    <- prep$cLLA
  ni      <- prep$ni

  # -----------------------------
  # Validate Bs
  # -----------------------------
  if (length(Bs) != n_LLA) {
    stop("Bs must be a list of length equal to n_LLA (number of longitudinal processes).")
  }

  # -----------------------------
  # Build surv.init using InitVal_LLAJEL_tv
  # -----------------------------
  surv.init <- InitVal_LLAJEL_tv(
    bBLUP      = bLLA,
    n_subject  = n_subject,
    n_LME      = n_LLA,
    Indcs      = Indcs,
    start      = start,
    stop       = stop,
    event      = event,
    W          = Wtime,
    ncw        = ncw,
    Wtime2     = Wtime2,
    ID_surv    = ID_surv,
    Bs         = Bs
  )

  B_list <- surv.init$B_list
  eta_n  <- surv.init$eta_n

  # -----------------------------
  # Run EM (LLA Time-Varying)
  # -----------------------------
  theta.refinedfastEM <- RefinedfastEM_LLA_tv(
    data         = train_dataset$LMM_dat,
    n_LLA        = n_LLA,
    ID_surv      = ID_surv,
    ni_surv      = ni_surv,
    Z.st         = Z.st,
    Y.st         = Y.st,
    Wker.st      = Wker.st,
    Wtime        = Wtime,
    Ysigma2      = Ysigma2,
    Bsigma       = Bsigma,
    bLLA         = bLLA,
    cLLA         = cLLA,
    n_subject    = n_subject,
    surv.init    = surv.init,
    ni           = ni,
    ncw          = ncw,
    B_list       = B_list,
    eta_n        = eta_n,
    gh.nodes     = gh.nodes,
    collect.hist = collect.hist,
    max.iter     = max.iter,
    tol          = tol,
    diff.type    = diff.type,
    post.process = post.process,
    verbose      = verbose
  )
  gc()

  converge <- if (!is.null(theta.refinedfastEM$history) && nrow(theta.refinedfastEM$history) < max.iter) 1 else 0

  # -----------------------------
  # Output (keep JEL-style structure)
  # -----------------------------
  result <- list()
  lamb.new <- as.data.frame(theta.refinedfastEM$coeffs$hazard)
  colnames(lamb.new) <- c("time", "bashaz")

  result$coefficients <- list(
    beta   = theta.refinedfastEM$coeffs$beta,
    phi    = theta.refinedfastEM$coeffs$phi,
    eta    = theta.refinedfastEM$coeffs$eta,
    Ysigma = sqrt(theta.refinedfastEM$coeffs$var.e),
    Bsigma = theta.refinedfastEM$coeffs$D,
    lamb   = lamb.new
  )

  result$logLik <- NULL
  result$call <- match.call()
  if (!is.null(theta.refinedfastEM$Vcov)) result$Vcov <- theta.refinedfastEM$Vcov
  result$est.bi <- theta.refinedfastEM$REs
  result$convergence <- if (converge == 1) "success" else "failure"
  if (!is.null(theta.refinedfastEM$postprocess.time)) result$time.SE <- theta.refinedfastEM$postprocess.time

  result$N <- sum(sapply(prep$Y, length))
  result$n <- n_subject
  result$d <- d
  result$dataMat <- list(B = theta.refinedfastEM$dmats$BB, ID_surv = ID_surv, Indcs = Indcs)

  class(result) <- "JEL"
  result
}