rm(list = ls())
library(MASS)
library(kernlab)
library(mvtnorm)
library(Matrix)

source('code/RBF.R')

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
