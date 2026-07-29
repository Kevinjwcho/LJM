
#========== Differentiate the S function with Forward Difference for Model I & II ==========#

PFDS_JEL <- function (model, theta, cvals, ncx, ncz, ncw, p, varNames, n, Z.st, Y.st, X.st, Ztime, nk, Wtime, Wtime2, Xtime, Xtime2, GQ,  rho, d, wGQ, ncz2, b, Ztime2.st,  X, Y, ID, N, Indcs, Z, B) {
  
  delta = cvals$delta;
  tol = min(cvals$tol.P, cvals$delta)/100;
  iter = cvals$max.iter 
  eta.name = varNames$eta.name;
  phi.names = varNames$phi.names
  beta.names = varNames$beta.names 
  
  n_eta = length(theta$eta)
  
  Index = Indcs$Index  
  Index0 = Indcs$Index0
  Index1 = Indcs$Index1
  Index2 = Indcs$Index2
  
  S <- Sfunc_JEL(model, theta, n, Z.st, Y.st, X.st, Ztime, nk, Wtime,  Wtime2, Xtime, Xtime2, GQ, Index, Index1, rho, d, wGQ, ncx, ncw, p, ncz, ncz2, b, Ztime2.st, Index0, X, Y, ID, N, Index2, Z, B)
  
  para <- List2Vec_JEL(theta)
  lamb.init <- theta$lamb
  len <- length(para)
  
  #===== Calculate the derivative of the score vector using forward difference =====#
  DS <- matrix(0, len, len)
  
  for (i in 1:len) {
    para1 <- para
    para1[i] <- para[i] + delta
    result <- LambGeneric_JEL(para1, lamb.init, tol, iter, ncz, ncx, ncw, n_eta, n, Z.st, Y.st, X.st, b, Ztime, Ztime2.st, nk, Wtime, Xtime, Wtime2, Xtime2, rho, Index0, Index1, Index, wGQ, model, GQ, d, Index2, B)
    para1.list <- Vec2List_JEL(para1, ncx, ncz, ncw, n_eta)
    theta.input1 <- list(beta = para1.list$beta, phi = para1.list$phi, eta = para1.list$eta, 
                         Ysigma = para1.list$Ysigma, Bsigma = para1.list$Bsigma, lamb = result$lamb)
    S1 <- Sfunc_JEL(model, theta.input1, n, Z.st, Y.st, X.st, Ztime, nk, Wtime, Wtime2, Xtime, Xtime2, GQ, Index, Index1, rho, d, wGQ, ncx, ncw, p, ncz, ncz2, b, Ztime2.st, Index0, X, Y, ID, N, Index2, Z, B)
    DS[i, ] <- (S1 - S) / delta
  }
  
  #========== make the DS matrix symmetric ==========#
  DS <- (DS + t(DS)) / 2 
  V <- -solve(DS);
  Veta.name <- paste("eta:", eta.name, sep = "")
  Vnames <- c(paste(rep("beta:", ncx), beta.names, sep = ""), paste(rep("phi:", ncw), phi.names, sep = ""),
              Veta.name, "sigma.e", paste("Bsigma.", 1:p, sep = ""))
  dimnames(V) <- list(Vnames, Vnames)
  
  return(V)
}
