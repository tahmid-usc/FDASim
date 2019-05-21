# RBF kernel

Jmat <- function(m) {
  return(matrix(rep(1, m^2), ncol = m))
}

BDmat <- function(repnum) {
  mat <- lapply(repnum, Jmat)
  as.matrix(bdiag(mat))
}

covmat <- function(trainx, repnum, theta) {
  rbf1 <- rbfdot(sigma = 1/theta[1])
  rbf2 <- rbfdot(sigma = 1/theta[3])
  k1 <- theta[2] * kernelMatrix(rbf1, x = trainx)
  k2 <- theta[4] * kernelMatrix(rbf2, x = trainx)
  k2 <- k2 * BDmat(repnum)
  return(k1 + k2)
}

testcov <- function(x, y, theta) {
  rbf1 <- rbfdot(sigma = 1/theta[1])
  #rbf2 <- rbfdot(sigma = abs(theta[3]))
  k1 <- theta[2] * kernelMatrix(rbf1, x = x, y = y)
  #k2 <- abs(theta[4]) * kernelMatrix(rbf2, x = x, y = y)
  return(k1)
}

testmat <- function(x, theta) {
  rbf1 <- rbfdot(sigma = 1/theta[1])
  rbf2 <- rbfdot(sigma = 1/theta[3])
  k1 <- theta[2] * kernelMatrix(rbf1, x = x)
  k2 <- theta[4] * kernelMatrix(rbf2, x = x)
  return(k1 + k2)
}