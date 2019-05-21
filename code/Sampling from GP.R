rm(list = ls())
library(MASS)
library(kernlab)
library(mvtnorm)
#library(Matrix)

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



# Generate mean function from GP
l1 <- 0.01
l2 <- 0.01

x <- 0:100
mu <- rep(0,101)
sigma.mu <- ker(x=x,l=l1, sigf=1)
mut <- mvrnorm(1,mu,sigma.mu)

# Generate random effect
sigma.gt <- ker(x=x,l=l2, sigf=1)
gt <- mvrnorm(1,mu,sigma.gt)

# Generate observed curve

y <- mut + gt #without measurement error
#y <- mut + gt + rnorm(101, 0 , .5) #with measurement error

plot(x, mut, lwd = 3, col = 1, type = 'l', ylab = 'Y')
lines(y, col='grey', lwd=3)

#---Plot non sparse

l1 <- 0.01
l2 <- 0.01

x <- 1:100
mu <- rep(0,100)
sigma.mu <- ker(x=x,l=l1, sigf=1)
sigma.gt <- ker(x=x,l=l2, sigf=1)
mut <- mvrnorm(1,mu,sigma.mu)

plot(x, mut, lwd = 3, col = 1, type = 'l', ylab = 'Y', ylim = c(-4,4))

for(i in 1:5) {
  gt <- mvrnorm(1,mu,sigma.gt)
  #y <- mut + gt 
  y <- mut + gt + rnorm(100, 0 , .08)
  lines(y, col= 'grey', lwd = 2)
}

legend('bottomright', c('Mean fucntion', 'Observations'),lty = 1, col = c(1, 'grey'), bty = 'n')


#---Plot sparse


plot(x, mut, lwd = 3, col = 1, type = 'l', ylab = 'Y', ylim = c(-4,4))

for(i in 1:5) {
  sel <- sort(sample(1:100, 10))
  gt <- mvrnorm(1,mu,sigma.gt)
  #y <- mut + gt 
  y <- mut + gt + rnorm(100, 0 , .08)
  lines(sel, y[sel], col= 'grey', lwd = 2, type = 'b')
}

legend('bottomright', c('Mean fucntion', 'Observations'),lty = 1, col = c(1, 'grey'), bty = 'n')

