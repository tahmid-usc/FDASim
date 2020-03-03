rm(list = ls())
library(MASS)
library(kernlab)
library(mvtnorm)
library(Matrix)
library(optimx)
library(fdapace)
library(dplyr)

source('code/RBF.R')
source('code/functionsNew.R')
#source('code/laplace.R')
#source('code/matern52.R')
#source('code/matern32.R')
#source('code/functions.R')
#source('code/multistart.R')


#-------------------------------------------------


# create training dataset

createfdt <- function(n) {
  muf <- function(x) {
    1 * sin(x * 10) + 1 * x
    #2 +  x
  }
  
  
  train0 <- fdagen(n = 50, maxt = 10, muf = muf, theta = c(.1,.1,.01))
  train0$z <- rep(1, dim(train0)[1])
  
  
  
  
  muf <- function(x) {
    #1 * sin(x * 10) + 1 * x
    x
  }
  
  train1 <- fdagen(n = 50, maxt = 10, muf = muf, theta = c(.1,.1,.01))
  train1$id <- train1$id + n
  train1$z <- rep(0, dim(train1)[1])
  
  
  
  train <- rbind(train0, train1)
  
  return(train)
}


train <- createfdt(50)
test <- createfdt(50)


#--------------------------------

train <- split(train,train$z)
fet <- lapply(train, feature)


#------------------


log.prob <- c()

for(i in unique(test$id)) {
  
  fit1 <- fit.gp(fet$`0`, testx = test[test$id == i, 1], testy = test[test$id == i, 2])
  fit2 <- fit.gp(fet$`1`, testx = test[test$id == i, 1], testy = test[test$id == i, 2])
  log.prob <- rbind(log.prob, c(fit1, fit2))
  
}

y <- c()
for (i in unique(test$id)) {
  y <- rbind(y, as.character((test[test$id == i, 4])[1]))
}

y.pred <- log.prob[,1] > log.prob[,2]
pred <- ifelse(y.pred == 1, 0, 1)
table(pred, y)
err <- mean(pred != y)
err


# plot


curve(sin(x * 10) + x, 0 , 1, lwd = 3)
curve(1 * x, 0 , 1, add = T, lwd = 3, col = 2)
for(i in 1:50) {
  lines(test[test$id == i, 1], test[test$id == i, 2], type = 'b', col = (test[test$id == i, 4] ) )
}
