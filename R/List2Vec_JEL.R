
List2Vec_JEL <- function (theta) {

  Bsigma <- theta$Bsigma
  Bsigma <- if (is.matrix(Bsigma)) Bsigma[lower.tri(Bsigma, diag = TRUE)] else Bsigma
  
  para <- c(theta$beta, theta$phi, theta$eta, theta$Ysigma, Bsigma)
  return(para)
}
