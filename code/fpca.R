rm(list = ls())
library(fdapace)
library(MASS)
library(kernlab)
library(mvtnorm)
library(Matrix)
library(optimx)

source('code/RBF.R')
source('code/laplace.R')
source('code/matern52.R')
source('code/matern32.R')
source('code/functions.R')
#source('code/multistart.R')

# Generate training data----

train <- gen(n = 50, theta1 = c(1000,1,.1), theta2 = c(1000,1,.1), maxt = 5)
train <- split(train, train$z)

# standardize

my <- mean(c(train$train1$y, train$train2$y))
stdy <- sd(c(train$train1$y, train$train2$y))

train$train1$y <- (train$train1$y - my) / stdy
train$train2$y <- (train$train2$y - my) / stdy

# Fit FPCA

dt1 <- MakeFPCAInputs(IDs = train$train1$id, train$train1$x, train$train1$y, sort = FALSE)
dt2 <- MakeFPCAInputs(IDs = train$train2$id, train$train2$x, train$train2$y, sort = FALSE)


pca1 <- FPCA(dt1$Ly, dt1$Lt, optns = list(dataType = 'Sparse', usergrid= T))
pca2 <- FPCA(dt2$Ly, dt2$Lt, optns = list(dataType = 'Sparse'))

# Fit GP

fet <- lapply(train, feature)

# Compare

a1 <- 2
b1 <- .1
a2 <- 2.5
b2 <- 1
c2 <- 0.06

fx1 <- function(t) { 
  a1 <- 2
  b1 <- .1
  return(((a1 + b1 * t)-my)/stdy) 
}

fx2 <- function(t) { 
  a2 <- 2.5
  b2 <- 1
  c2 <- 0.06
  return(((a2 + b2 * sin(t/5) + c2 * t)-my)/stdy) 
}

mu1 <- (mu1 - my) / stdy
mu2 <- a2 + b2 * sin(pca2$workGrid/5) + c2 * pca2$workGrid 
mu2 <- (mu2 - my) / stdy

mu1.gp <- gpsmooth(pca1$workGrid, fet$train1)
mu2.gp <- gpsmooth(pca2$workGrid, fet$train2)

#FPCA
mean((mu1 - pca1$mu)^2)
mean((mu2 - pca2$mu)^2)
# GP
mean((mu1 - mu1.gp)^2)
mean((mu2 - mu2.gp)^2)

# MISE ---------------------------

residfunc <- function(t, fp, gp, fx) {
  ferr <- (fpcamu(fp, t) - fx(t))^2
  gperr <- (gpsmooth(t, gp) - fx(t))^2
  return(list(fp = ferr, gp = gperr))
}

xt <- 0:100
residfunc(xt, fp = pca1, gp = fet$train1, fx = fx1)
baal <- function(t) {residfunc(t, fp = pca1, gp = fet$train1, fx = fx1)$fp}
integrate(baal, 0 , 100, subdivisions = 10000)


#-------------------------

plot(pca1$workGrid, mu1, type = 'l', lwd = 4, main = 'Mean function estimation (Linear)',
     xlab = 'x', ylab = 'y')
lines(pca1$workGrid, pca1$mu, lty = 2, lwd = 3, col = 2)
lines(pca1$workGrid, mu1.gp, lty = 2, lwd = 3, col = 3)
legend('bottomright', c('True', 'FPCA', 'GP'), col = 1:3, lty = 1, bty = 'n')


plot(pca2$workGrid, mu2, type = 'l', lwd = 4, main = 'Mean function estimation (Periodic)',
     xlab = 'x', ylab = 'y')
lines(pca2$workGrid, pca2$mu, lty = 2, lwd = 3, col = 2)
lines(pca2$workGrid, mu2.gp, lty = 2, lwd = 3, col = 3)
legend('bottomright', c('True', 'FPCA', 'GP'), col = 1:3, lty = 1, bty = 'n')

