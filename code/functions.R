
# Generate data fromm linear population
gen <- function(n = 100, maxt = 20, theta1 = rep(1,3), theta2 = rep(1,3)) {
  a1 <- 2
  b1 <- .1
  a2 <- 2.5
  b2 <- 1
  c2 <- 0.06
  train1 <- data.frame()
  train2 <- data.frame()
  #n number of functions in the sample data
  n.time <- sample(1:maxt, size=n, replace=T) 
  
  for(i in 1:n) {
    id <- rep(i, n.time[i])
    x <- sort(runif(n.time[i], 0 , 100))
    mu <- a1 + b1 * x
    gt <- mvrnorm(1, rep(0, n.time[i]), ker(x, l = theta1[1], sigf = theta1[2])) 
    y <- mu + gt + rnorm(n.time[i], 0, theta1[3])
    train1 <- rbind(train1, cbind(x, y, id))
  }
  
  z <- rep('train1', dim(train1)[1])
  train1 <- cbind(train1,z)
  
  n.time <- sample(1:maxt, size=n, replace=T) 
  
  for(i in 1:n) {
    id <- rep(i, n.time[i]) + n
    x <- sort(runif(n.time[i], 0 , 100))
    mu <- a2 + b2 * sin(x/5) + c2 * x
    gt <- mvrnorm(1, rep(0, n.time[i]), ker(x, l = theta2[1], sigf = theta2[2])) 
    y <- mu  + gt + rnorm(n.time[i], 0, theta2[3])
    train2 <- rbind(train2, cbind(x, y, id))
  }
  z <- rep('train2', dim(train2)[1])
  train2 <- cbind(train2,z)
  
  train <- rbind(train1, train2)
  return(train)
}


# Extract relevant feature of the data

feature <- function(train) {
  
  n <- length(unique(train$id))
  N <- dim(train)[1]
  trainx <- train$x
  trainy <- train$y
  rep <- as.numeric(table(train$id))
  hyper <- Hyper(trainx = trainx, trainy = trainy, repnum = rep, N = N)
  kinv <- chol2inv(chol(covmat(trainx, rep , theta = hyper[1:4]) + hyper[5] * diag(N)))
  return(list(n = n,N = N,rep = rep, trainx = trainx, trainy = trainy, hyper = hyper, kinv = kinv))
  
}


# Estimate hypoerparameter

Hyper <- function(trainx, trainy, repnum, N) {
  
  marlik <- function(theta) {
    theta <- theta^2
    kxx <- covmat(trainx = trainx, repnum = repnum, theta = theta[1:4])
    -dmvnorm(x = trainy, sigma = kxx + theta[5] * diag(N), log = T)
  }
  
  hyp <- optim(par=rep(1, 5), fn = marlik, method = 'Nelder-Mead',
               control=list(maxit = 100000))
  print(hyp)
  return(hyp$par^2)
  
}


# Fit a smooth curve--------------------------------

gpsmooth <- function(x, trainlist) {
  
  kxx <- covmat(trainx = trainlist$trainx, repnum = trainlist$rep, 
                theta = trainlist$hyper[1:4])
  kx <-  testcov(x = x, y = trainlist$trainx, theta = trainlist$hyper[1:4])
  kinv <-  chol2inv(chol(kxx + trainlist$hyper[5] * diag(trainlist$N)))
  k <- kx %*% kinv
  pred <- k %*% as.matrix(trainlist$trainy)
  return(pred)
  
}


# Posterior mean and cov matrix for test curve

fit.gp <- function(train, testx, testy) {
  
  n <- length(testx)
  #kxx <- covmat(trainx, repnum , theta = theta[1:4])
  kx <-  testcov(x = testx, y = train$trainx, theta = train$hyper[1:4])
  k <- kx %*% train$kinv
  mu <- k %*% as.matrix(train$trainy)    
  sigma <- testmat(x = testx, theta = train$hyper[1:4]) + (train$hyper[5] * diag(n)) - (k %*% t(kx))
  sigma <- as.matrix(forceSymmetric(sigma))
  logprob <- dmvnorm(x = testy, mean = mu, sigma = sigma, log = T)
  return(logprob)
}



# Local linear smoother (estimate of FPCA mean)
fpcamu <- function(fpca, t) {
  
  to <- order(t)
  t <- sort(t)
  xt <- unlist(fpca$inputData$Lt)
  yt <- unlist(fpca$inputData$Ly)
  dt <- data.frame(xt, yt)
  dt <- dt[order(dt$xt),]
  
  mu <- Lwls1D(bw = fpca$bwMu, kernel_type = 'gauss', xin = dt$xt, 
               yin = dt$yt, xout = t)
  mu <- mu[to]
  return(mu)
}
