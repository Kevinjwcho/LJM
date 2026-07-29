# Other functions ---------------------------------------------------------
vech <- function(x) x[lower.tri(x, diag = T)]
tr <- function(x) sum(diag(x))
an <- as.numeric
repCols <- function(X, n = 3) X[,rep(1:2, n)] 
repVec <- function(x, n = 3) rep(x, n)
