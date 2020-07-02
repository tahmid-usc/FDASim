# Matern 5/2

ker <- function(x, l = 1, sigf = 1) {
  r <- as.matrix(dist(x))
  k <- sigf * (1 + (sqrt(5) * r / l) + ((5 * r^2) / (3 * l^2))) * exp(- sqrt(5) * r / l)
  return(k)
}

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


Jmat <- function(m) {
  return(matrix(rep(1, m^2), ncol = m))
}

BDmat <- function(repnum) {
  mat <- lapply(repnum, Jmat)
  as.matrix(bdiag(mat))
}

covmat <- function(trainx, repnum, theta) {
  
  k1 <- covf(x = trainx, l = theta[1], sigf = theta[2])
  #k2 <- covf(x = trainx, l = theta[3], sigf = theta[4])
  rbf2 <- rbfdot(sigma = 1/theta[3])
  k2 <- theta[4] * kernelMatrix(rbf2, x = trainx)
  k2 <- k2 * BDmat(repnum)
  return(k1 + k2)
}

testcov <- function(x, y, theta) {
  
  k1 <- covfTest(x = x, y = y, l = theta[1], sigf = theta[2])
  return(k1)
}

testmat <- function(x, theta) {
  
  k1 <- covf(x, l = theta[1], sigf = theta[2])
  #k2 <- covf(x, l = theta[3], sigf = theta[4])
  rbf2 <- rbfdot(sigma = 1/theta[3])
  k2 <- theta[4] * kernelMatrix(rbf2, x = trainx)
  return(k1 + k2)
}