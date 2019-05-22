
# Generate data fromm linear population
gen1 <- function(n = 100, maxt = 20, theta = rep(1,3)) {
  a <- 2
  b <- .1
  #n number of functions in the sample data
  n.time <- sample(1:maxt, size=n, replace=T) 
  train <- data.frame()
  
  for(i in 1:n) {
    id <- rep(i, n.time[i])
    x <- sort(runif(n.time[i], 0 , 100))
    mu <- a + b * x
    gt <- mvrnorm(1, rep(0, n.time[i]), ker(x, l = theta[1], sigf = theta[2])) 
    y <- mu + gt #+ rnorm(n.time[i], 0, theta[3])
    train <- rbind(train, cbind(x, y, id))
  }
  z <- rep(1, dim(train)[1])
  train <- cbind(train,z)
  return(train)
}


# Generate data from periodic population

gen2 <- function(n = 100, maxt = 20, theta = rep(1,3)) {
  a <- 2.5
  b <- 1
  c <- 0.06
  #n number of functions in the sample data
  n.time <- sample(1:maxt, size=n, replace=T) 
  train <- data.frame()
  
  for(i in 1:n) {
    id <- rep(i, n.time[i])
    x <- sort(runif(n.time[i], 0 , 100))
    mu <- a + b * sin(x/5) + c * x
    gt <- mvrnorm(1, rep(0, n.time[i]), ker(x, l = theta[1], sigf = theta[2])) 
    y <- mu  + gt + rnorm(n.time[i], 0, theta[3])
    train <- rbind(train, cbind(x, y, id))
  }
  z <- rep(2, dim(train)[1])
  train <- cbind(train,z)
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
  
  hyp <- optim(par=rep(1, 5), fn = marlik, method = 'BFGS',
               control=list(maxit = 10000, ndeps = rep(.00001, 5)))
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

