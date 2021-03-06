
# Generate data 

# Mean function 
muf1 <- function(x) {
  fx <- #.25 * dnorm(x, mean = .2, sd = .05)
    dnorm(x, mean = .2, sd = .05) + 
    dnorm(x, mean = .5, sd = .03) +   
    dnorm(x, mean = .8,sd = .02) 
  return(fx)
}


muf2 <- function(x) {
  fx <- 
     ddoublex(x, .2, .05) + 
     ddoublex(x, .5,  .03) +   
     ddoublex(x, .8, .02) 
  return(fx)
}


# data frame

fdagen <- function(n = 100, maxt = 20, muf, theta = rep(1,3)) {
  
  train <- data.frame()
  #n number of functions in the sample data
  n.time <- sample(2:maxt, size=n, replace=T) 
  
  for(i in 1:n) {
    id <- rep(i, n.time[i])
    x <- sort(runif(n.time[i], 0 , 1))
    mu <- muf(x)
    gt <- mvrnorm(1, rep(0, n.time[i]), ker(x, l = theta[1], sigf = theta[2])) 
    #y <- mu + gt
    y <- mu + gt + rnorm(n.time[i], 0, theta[3])
    train <- rbind(train, cbind(x, y, id))
  }
  return(train)
}


# Extract relevant feature of the data

feature <- function(train) {
  
  n <- length(unique(train$id))
  N <- dim(train)[1]
  trainx <- train$x
  trainy <- train$y
  # standardize training data
  ymean <- mean(trainy)
  ysd <- sd(trainy)
  trainy <- (trainy - ymean)/ysd
  
  rep <- as.numeric(table(train$id))
  hyper <- Hyper(trainx = trainx, trainy = trainy, repnum = rep, N = N)
  kinv <- chol2inv(chol(covmat(trainx, rep , theta = hyper[1:4]) + hyper[5] * diag(N)))
  return(list(n = n,N = N,rep = rep, trainx = trainx, trainy = trainy, hyper = hyper, kinv = kinv, ymean = ymean, ysd = ysd))

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
  #kinv <-  chol2inv(chol(kxx + trainlist$hyper[5] * diag(trainlist$N)))
  k <- kx %*% trainlist$kinv
  pred <- k %*% as.matrix(trainlist$trainy)
  pred <-  trainlist$ymean + trainlist$ysd * pred
  return(pred)
  
}



credible_band <- function(x, trainlist) {
  
  kxx <- covmat(trainx = trainlist$trainx, repnum = trainlist$rep, 
                theta = trainlist$hyper[1:4])
  kx <-  testcov(x = x, y = trainlist$trainx, theta = trainlist$hyper[1:4])
  #kinv <-  chol2inv(chol(kxx + trainlist$hyper[5] * diag(trainlist$N)))
  k <- kx %*% trainlist$kinv
  pred <- k %*% as.matrix(trainlist$trainy)
  pred <-  trainlist$ymean + trainlist$ysd * pred
  
  n <- length(x)
  sigma <- testmat(x = x, theta = trainlist$hyper[1:4]) + (trainlist$hyper[5] * diag(n)) - (k %*% t(kx))
  sigma <- as.matrix(forceSymmetric(sigma))
  
  return(list(pred = pred, sigma = sigma))
  
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
fpcamu <- function(t, fpca) {
  
  to <- order(t)
  t <- sort(t)
  xt <- unlist(fpca$inputData$Lt)
  yt <- unlist(fpca$inputData$Ly)
  dt <- data.frame(xt, yt)
  dt <- arrange(dt, xt)
  
  mu <- Lwls1D(bw = fpca$bwMu, kernel_type = 'gauss', xin = dt$xt, 
               yin = dt$yt, xout = t)
  mu <- mu[to]
  return(mu)
}




## MISE ---------------------------

residfunc <- function(x, muf, est, est_arg) {
  ((muf(x) - est(x, est_arg))^2)
}


# generate data


funcgen <- function(muf, theta) {
  x <- seq(0,1,length.out = 100)
  mu <- muf(x)
  y <- mu + mvrnorm(1, rep(0,100), ker(x, l = theta[1], sigf = theta[2])) 
  return(list(x=x, y = y))
}


# RMSE

ase <- function(n, method, muf, est, est_arg) {
  if(method == 'mc') time <- runif(n, 0, 1) else time <- seq(0, 1, length.out = n)
  ase <-  mean((residfunc(time, muf, est, est_arg)))
  return(ase)
}

#vase <- Vectorize(ase, vectorize.args = 'n')
