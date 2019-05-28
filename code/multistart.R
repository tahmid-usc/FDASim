

Hyper <- function(trainx, trainy, repnum, N) {
  
  marlik <- function(theta) {
    theta <- theta^2
    kxx <- covmat(trainx = trainx, repnum = repnum, theta = theta[1:4])
    -dmvnorm(x = trainy, sigma = kxx + theta[5] * diag(N), log = T)
  }
  
  parval <- c(.01, .05, .1, .5, 1, 2, 3, 4, 5)
  parmat <- matrix(parval, nrow = length(parval), ncol = 5)
  hyp <- multistart(parmat, fn = marlik, method='BFGS', control = list(maxit = 10000))
  hyp <- as.numeric(hyp[which(hyp$value == min(hyp$value)), 1:5])^2
  return(hyp)
  
}
