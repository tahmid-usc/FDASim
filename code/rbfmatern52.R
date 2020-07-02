# RBF kernel

ker <- function(x, l, sigf) {
  rbf <- rbfdot(sigma = 1/l)
  return(sigf * kernelMatrix(rbf, x = x))
}

Jmat <- function(m) {
  return(matrix(rep(1, m^2), ncol = m))
}

BDmat <- function(repnum) {
  mat <- lapply(repnum, Jmat)
  as.matrix(bdiag(mat))
}

#---


covf <- function(x, l = 1, sigf = 1) {
  r <- as.matrix(dist(x))
  k <- sigf * (1 + (sqrt(5) * r / l) + ((5 * r^2) / (3 * l^2))) * exp(- sqrt(5) * r / l)
  return(k)
}


covfTest <- function(x, y, l = 1, sigf = 1) {
  r <- abs(outer(x, y, "-"))
  k <- sigf * (1 + (sqrt(5) * r / l) + ((5 * r^2) / (3 * l^2))) * exp(- sqrt(5) * r / l)
  return(k)
}

#---


covmat <- function(trainx, repnum, theta) {
  rbf1 <- rbfdot(sigma = 1/theta[1])
  k1 <- theta[2] * kernelMatrix(rbf1, x = trainx)
  k2 <- covf(x = trainx, l = theta[3], sigf = theta[4])
  k2 <- k2 * BDmat(repnum)
  return(k1 + k2)
}

testcov <- function(x, y, theta) {
  rbf1 <- rbfdot(sigma = 1/theta[1])
  k1 <- theta[2] * kernelMatrix(rbf1, x = x, y = y)
  return(k1)
}

testmat <- function(x, theta) {
  rbf1 <- rbfdot(sigma = 1/theta[1])
  k1 <- theta[2] * kernelMatrix(rbf1, x = x)
  k2 <- covf(x = trainx, l = theta[3], sigf = theta[4])
  return(k1 + k2)
}