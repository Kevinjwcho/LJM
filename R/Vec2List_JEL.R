Vec2List_JEL <- function (para, ncx, ncz, ncw, n_eta) {

  p <- ncz * (ncz + 1) / 2
  beta <- para[1 : ncx]
  phi <- if (ncw > 0) para[(ncx + 1) : (ncx + ncw)] else numeric(0)
  eta <- para[(ncx + ncw + 1):(ncx + ncw + n_eta)]
  Ysigma <- para[ncx + ncw + n_eta + 1]
  Bsigma <- para[(ncx + ncw + n_eta + 2):(ncx + ncw + p + 1+ n_eta)]
  if (ncz > 1) {
    Bsigma.new <- matrix(0, ncz, ncz)
    Bsigma.new[lower.tri(Bsigma.new, diag = TRUE)] <- Bsigma
    Bsigma.new <- Bsigma.new + t(Bsigma.new) - diag(diag(Bsigma.new))
    Bsigma <- Bsigma.new
  }
  result <- list(beta = beta, phi = phi, eta = eta, Ysigma = Ysigma, Bsigma = Bsigma)
  return(result)
}

