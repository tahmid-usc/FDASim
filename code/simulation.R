rm(list = ls())
library(MASS)
library(kernlab)
library(mvtnorm)
library(Matrix)

#Kernel------------------------

ker <- function (x, l, sigf) {
  rbf <- rbfdot(sigma = l)
  sigf * kernelMatrix(rbf, x)
} 

#Laplace 

ker <- function (x, l, sigf) {
  rbf <- laplacedot(sigma = l)
  sigf * kernelMatrix(rbf, x)
} 


# Generate data with/without measurement error with linear trend-------

a1 <- 2
b1 <- .1

x <- 0:100
mu <- rep(0,101)
sigma <- ker(x=x,l=.011, sigf=1)

ft.lin <- a1 + b1 * x
gt <- mvrnorm(1,mu,sigma)
#y.lin <- ft.lin + gt
y.lin <- ft.lin + gt + rnorm(101, 0 , .5)

curve(a1 + b1 * x, 0 , 100, lwd = 3, col = 'grey')
lines(y.lin, col=1, lwd=3)


# Generate data with/without measurement error with periodic trend-------


a2 <- 2
b2 <- 1
c2 <- .1

x <- 0:100
mu <- rep(0,101)
sigma <- ker(x=x,l=.012, sigf=1)

ft.per <-  a2 + b2 * sin(x/5) + c2 * x
gt <- mvrnorm(1,mu,sigma)
y.per <- ft.per + gt
#y.per <- ft.per + gt + rnorm(101, 0 , .5)

curve(a2 + b2 * sin(x/5) + c2 * x, 0 , 100, lwd = 3, col = 1)
lines(y.per, col= 'grey', lwd=3)


#---Plot

a1 <- 2
b1 <- .1

x <- 0:100
mu <- rep(0,101)
sigma <- ker(x=x,l=.005, sigf=1)

curve(a1 + b1 * x, 0, 100, lwd = 4, ylab = 'Y', main = 'Sample functions from linear 
      population (Gaussian Kernel) with measurement error')

for(i in 1:10) {
  gt <- mvrnorm(1,mu,sigma)
  #y.lin <- ft.lin + gt
  y.lin <- ft.lin + gt + rnorm(101, 0 , .5)
  lines(y.lin, col= 'grey', lwd = 3)
}

#------------

a2 <- 2
b2 <- 1
c2 <- .1

x <- 0:100
mu <- rep(0,101)
sigma <- ker(x=x,l=.012, sigf=1)

ft.per <-  a2 + b2 * sin(x/5) + c2 * x

curve(a2 + b2 * sin(x/5) + c2 * x, 0, 100, lwd = 4, ylab = 'Y', main = 'Sample functions from 
periodic population (Gaussian Kernel) with measurement error')

for(i in 1:10) {
  gt <- mvrnorm(1,mu,sigma)
  #y.per <- ft.per + gt
  y.per <- ft.per + gt + rnorm(101, 0 , .5)
  lines(y.per, col= 'grey', lwd = 3)
}


#Generate training data----

# Linear mean
n.sam <- 100 #number of functions in the sample data
max.time <- 20
n.time.x1 <- sample(1:max.time, size=n.sam, replace=T) 
X1 <- list()
Y1 <- list()

a1 <- 2
b1 <- .1

for(i in 1:n.sam) {
  x <- sort(runif(n.time.x1[i], 0 , 100))
  mu <- a1 + b1 * x
  y <- mu + mvrnorm(1,rep(0,n.time.x1[i]),ker(x, l=.001, sigf=1)) 
  X1[[i]] <- x 
  Y1[[i]] <- y + rnorm(n.time.x1[i], 0, .1)
}

# Visualize ------
curve(a1 + b1 * x, 0 , 100, lwd=4, main = 'Linear Population')
for (i in 1:10) {
  #  lines(X1[[i]], Y1[[i]], lwd = .3, col = sample(1:30, 30, T))
  lines(X1[[i]], Y1[[i]], lwd = 2, col = 'grey', type = 'b', pch = 16)
  #browser()
}

# Perdiodic mean function ---------

X2 <- list()
Y2 <- list()

n.time.x2 <- sample(1:max.time, size=n.sam, replace=T) 

a2 <- 2.5
b2 <- 1
c2 <- .06

for(i in 1:n.sam) {
  x <- sort(runif(n.time.x2[i], 0 , 100)) 
  mu <- a2 + b2 * sin(x/5) + c2 * x
  y <- mu + mvrnorm(1, rep(0,n.time.x2[i]), ker(x, l=.001, sigf= .1)) 
  X2[[i]] <- x 
  Y2[[i]] <- y +  rnorm(n.time.x2[i], 0, .1)
}

# Visualize ------
curve(a2 + b2 * sin(x/5) + c2 * x, 0 , 100, lwd=4, main = 'Periodic Population')
for (i in 1:10) {
  lines(X2[[i]], Y2[[i]], lwd = 2, col = 'grey', type = 'b', pch = 16)
  #browser()
}

#Generate test data----

n.sam <- 1000 #number of functions in the sample data
max.time <- 20
n.time <- sample(1:max.time, size=n.sam, replace=T) 
X <- list()
Y <- list()

for(i in 1:(n.sam/2)) {
  x <- sort(runif(n.time[i], 0 , 100))
  mu <- a1 + b1 * x
  y <- mu + mvrnorm(1, rep(0, n.time[i]), ker(x, l=.001, sigf=1)) 
  X[[i]] <- x 
  Y[[i]] <- y + rnorm(n.time[i], 0, .1)
}


for(i in ((n.sam/2) + 1 ):n.sam) {
  x <- sort(runif(n.time[i], 0 , 100)) 
  mu <- a2 + b2 * sin(x/5) + c2 * x
  y <- mu + mvrnorm(1, rep(0, n.time[i]), ker(x, l=.001, sigf= .1)) 
  X[[i]] <- x 
  Y[[i]] <- y +  rnorm(n.time[i], 0, .1)
}

#Merge training data

train.x1 <- unlist(X1)
train.y1 <- unlist(Y1)

train.x2 <- unlist(X2)
train.y2 <- unlist(Y2)

#m <- mean(c(train.y1, train.y2))
m <- 0
train.y1 <- train.y1 - m
train.y2 <- train.y2 - m

# Estimating Hyperparameters

N.mal <- length(train.x1)
N.fem <- length(train.x2)

marlik.mal <- function(theta) {
  kxx.mal <- covmat(trainx = train.x1, repnum = n.time.x1, theta = theta[1:4])
  -dmvnorm(x = train.y1, 
           sigma = kxx.mal  + abs(theta[5]) * diag(N.mal), log = T)
}


hyp.mal <- optim(par=c(1,1,1,1,1), fn = marlik.mal, method = 'Nelder-Mead', 
                 control=list(maxit = 10000))
hyp.mal <- abs(hyp.mal$par)
hyp.mal


marlik.fem <- function(theta) {
  kxx.fem <- covmat(trainx = train.x2, repnum = n.time.x2, theta = theta[1:4])
  -dmvnorm(x = train.y2, 
           sigma = kxx.fem  + abs(theta[5]) * diag(N.fem), log = T)
}


hyp.fem <- optim(par=c(1,1,1,1,1), fn = marlik.fem, method = 'Nelder-Mead', 
                 control=list(maxit = 10000))

hyp.fem <- abs(hyp.fem$par)
hyp.fem

# Plotting predictions

kxx.mal <- covmat(trainx = train.x1, repnum = n.time.x1, theta = hyp.mal[1:4])
kx.mal <-  testcov(x = 0:100, y = train.x1, theta = hyp.mal[1:4])
kinv.mal  <-  chol2inv(chol(kxx.mal  + hyp.mal[5] * diag(N.mal)))
k.mal <- kx.mal %*% kinv.mal
pred.mal <- k.mal %*% as.matrix(train.y1)

kxx.fem <- covmat(trainx = train.x2, repnum = n.time.x2, theta = hyp.fem[1:4])
kx.fem <-  testcov(x = 0:100, y = train.x2, theta = hyp.fem[1:4])
kinv.fem  <-  chol2inv(chol(kxx.fem  + hyp.fem[5] * diag(N.fem)))
k.fem <- kx.fem %*% kinv.fem
pred.fem <- k.fem %*% as.matrix(train.y2)

curve(a1 + b1 * x - m, 0 , 100, lwd= 4, col = 1, lty = 2, 
      main = 'Posterior means for the pooled functiona data', ylab = 'Y')
lines(0:100, apply(matrix(0:100),2,function(x)(a2 + b2 * sin(x/5) + c2 * x - m)), lwd= 4, 
      col = 'grey', lty = 2)
lines(0:100, pred.mal, col = 1, lwd = 3)
points(train.x1, train.y1, pch = 16, cex= .7, col =1)
lines(0:100, pred.fem, col = 'grey', lwd = 3)
points(train.x2, train.y2, pch = 16, cex= .7, col = 'grey')
legend('bottomright', c('Class 1', 'Class 2', 'True Class 1', 'True Class 2'), 
       col = c(1,'grey',1,'grey'), lty=c(1,1,2,2), bty = 'n')


# Posterior means for test curves

fit.gp <- function(test.x) {
  n <- length(test.x)
  kx.mal <- testcov(x = test.x, y = train.x1, theta = hyp.mal[1:4])
  k.mal <- kx.mal %*% kinv.mal
  mu.mal <- k.mal %*% as.matrix(train.y1)
  sigma.mal <- testmat(x = test.x, theta = hyp.mal[1:4]) + diag(n) - (k.mal %*% t(kx.mal))
  sigma.mal <- as.matrix(forceSymmetric(sigma.mal))
  
  kx.fem <- testcov(x = test.x, y = train.x2, theta = hyp.fem[1:4])
  k.fem <- kx.fem %*% kinv.fem
  mu.fem <- k.fem %*% as.matrix(train.y2)
  sigma.fem <- testmat(x = test.x, theta = hyp.fem[1:4]) + diag(n) - (k.fem %*% t(kx.fem))
  sigma.fem <- as.matrix(forceSymmetric(sigma.fem))
  
  return(list(pred.mal = mu.mal, pred.fem = mu.fem, sigma.mal = sigma.mal, sigma.fem = sigma.fem))
}

i <- sample(1:n.sam, 1)
fit <- fit.gp(test.x = X[[i]])
loglik.mal <- dmvnorm((Y[[i]] - m), mean = fit$pred.mal, sigma = fit$sigma.mal, log = T)
loglik.fem <- dmvnorm((Y[[i]] - m), mean = fit$pred.fem, sigma = fit$sigma.fem, log = T)

color <- ifelse(i < ((n.sam/2)+1) , 1 , 2)
curve(a1 + b1 * x - m, 0 , 100, lwd= 3, col = 1, lty = 2, xlab = 'X', ylab = 'Y')
lines(0:100, apply(matrix(0:100),2,function(x)(a2 + b2 * sin(x/5) + c2 * x - m)), lwd= 3, col = 2, lty = 2)
lines(X[[i]], fit$pred.mal, col = 1, lwd = 3)
lines(X[[i]], fit$pred.fem, col = 2, lwd = 3)
lines(X[[i]], Y[[i]], col = 3, lwd = 2)
text(20,11, paste('Class', color, '\n', 'Class 1 loglik',  round(loglik.mal,4), '\n', 'Class 2 loglik', round(loglik.fem, 4)))
legend('bottomright', c('Class 1', 'Class 2', 'Test'), col = 1:3, lty=1, bty = 'n')

#---Classifying (USE this one)-------------

log.prob <- c()

for(i in 1:n.sam) {
  
  fit <- fit.gp(test.x = X[[i]])
  log.prob <- rbind(log.prob, 
                    c(dmvnorm((Y[[i]] - m), mean = fit$pred.mal, sigma = fit$sigma.mal, log = T), 
                      dmvnorm((Y[[i]] - m), mean = fit$pred.fem, sigma = fit$sigma.fem, log = T)))
}

#cbind(log.prob, log.prob[,1]>log.prob[,2], rep(c(1,0),each=(n.sam/2)) )
pred <- as.numeric(log.prob[,1]>log.prob[,2])
true <- rep(c(1,0), each = (n.sam/2))
table(pred, true)
mean(pred != true)


# See missclassifier ones
miss <- (1:n.sam)[pred != true]
cbind(log.prob[miss,],true[miss], pred[miss])

for (i in miss) {
  fit <- fit.gp(test.x = X[[i]])
  loglik.mal <- dmvnorm((Y[[i]] - m), mean = fit$pred.mal, sigma = fit$sigma.mal, log = T)
  loglik.fem <- dmvnorm((Y[[i]] - m), mean = fit$pred.fem, sigma = fit$sigma.fem, log = T)
  
  color <- ifelse(i < ((n.sam/2)+1) , 1 , 2)
  true <- ifelse(color == 1, 'Linear', 'Periodic')
  curve(a1 + b1 * x - m, 0 , 100, lwd= 3, col = 1, lty = 2, xlab = 'X', ylab = 'Y')
  lines(0:100, apply(matrix(0:100),2,function(x)(a2 + b2 * sin(x/5) + c2 * x - m)), lwd= 3, col = 2, lty = 2)
  lines(X[[i]], fit$pred.mal, col = 1, lwd = 3, type = 'b')
  lines(X[[i]], fit$pred.fem, col = 2, lwd = 3, type = 'b')
  lines(X[[i]], Y[[i]], col = 3, lwd = 2, type = 'b')
  text(20,11, paste( true, '\n', 'Linear loglik',  round(loglik.mal,4), '\n', 'Periodic loglik', round(loglik.fem, 4)))
  legend('bottomright', c('Linear', 'Periodic', 'Test'), col = 1:3, lty=1, bty = 'n')
  
} 





