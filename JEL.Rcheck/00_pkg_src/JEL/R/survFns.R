surv.mod <- function(cph, data, l0.init = NULL){
  uids <- unique(data$id)
  if("coxph.null"%in%class(cph)) message("Null model")
  # Survfit
  sf <- summary(survfit(cph))

  # initialise empty stores
  Fi <- matrix(NA, nr = length(uids), nc = 2)
  Di <- l0i <- numeric(length(uids))
  Fu <- l0u <- surv.times <- list()

  if(is.null(l0.init)){
    l0 <- diff(c(0, sf$cumhaz))
  }else{
    l0 <- l0.init
  }

  # loop
  for(i in 1:length(uids)){
    # ID = data$id
    ID = uids[i]
    i.dat <- subset(data, id == ID)
    Di[i] <- unique(i.dat$event)
    stop_i <- unique(i.dat$stop)
    surv.times[[i]] <- which(sf$time <= stop_i)
    Fi[i,] <- c(1, stop_i)
    Fu[[i]] <- cbind(1, rep(1, length(which(sf$time <= stop_i))))
    l0u[[i]] <- l0[which(sf$time <= stop_i)]
    # Robust: avoid floating-point `sf$time == stop_i` (was returning integer(0)
    # when init Cox is rank-deficient, triggering "replacement has length zero").
    # Use match() with NA-safe fallback, and accept a tiny tolerance.
    if (Di[i] == 1) {
      idx <- match(stop_i, sf$time)
      if (is.na(idx)) {
        diffs <- abs(sf$time - stop_i)
        if (length(diffs) > 0 && min(diffs) < 1e-8) {
          idx <- which.min(diffs)
        }
      }
      l0i[i] <- if (is.na(idx) || length(idx) == 0L) 0 else l0[idx]
    } else {
      l0i[i] <- 0
    }
    # Check if censored before first failure time
    if(Di[i] == 0 & unique(i.dat$stop) <= min(sf$time)){ l0u[[i]] <- 0; Fu[[i]] <- cbind(0, 0) }
  }

  nev <- c(); surv.ids <- list()
  p <- 1
  for(i in sf$time){
    nev[p] <- length(unique(data[which(data$stop == i),]$id))
    surv.ids[[p]] <- unique(data[which(data$stop >= i),]$id)
    p <- p+1
  }

  # output
  return(list(
    ft = sf$time,
    l0 = l0,
    nev = nev,
    surv.ids = surv.ids,
    surv.times = surv.times,
    l0i = l0i,
    Di = Di,
    Fi = Fi,
    Fu = Fu,
    l0u = l0u
  ))
}
