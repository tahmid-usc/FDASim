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



# Function create training dataset

createfdt <- function(n) {
  muf <- function(x) {
    1 * sin(x * 10) + 1 * x
  }
  
  
  train.per <- fdagen(n = 50, maxt = 10, muf = muf, theta = c(.1,.1,.08))
  train.per$z <- rep('Periodic', dim(train.per)[1])
  
  
  
  
  muf <- function(x) {
    0 + 1.2 * x
  }
  
  train.lin <- fdagen(n = 50, maxt = 10, muf = muf, theta = c(.1,.1,.08))
  train.lin$id <- train.lin$id + n
  train.lin$z <- rep('Linear', dim(train.lin)[1])
  
  
  
  train <- rbind(train.per, train.lin)
  
  return(train)
}



missclassError <- function(n1 = 50, n2 = 50) {
  train <- createfdt(n1)
  test <- createfdt(n2)
  
  
  #--------------------------------
  
  train <- split(train,train$z)
  fet <- lapply(train, feature)
  
  
  #------------------
  
  
  log.prob <- c()
  
  for(i in unique(test$id)) {
    
    fit1 <- fit.gp(fet$Periodic, testx = test[test$id == i, 1], testy = test[test$id == i, 2])
    fit2 <- fit.gp(fet$Linear, testx = test[test$id == i, 1], testy = test[test$id == i, 2])
    log.prob <- rbind(log.prob, c(fit1, fit2))
    
  }
  
  y <- c()
  for (i in unique(test$id)) {
    y <- rbind(y, as.character((test[test$id == i, 4])[1]))
  }
  
  y.pred <- log.prob[,1] > log.prob[,2]
  pred <- ifelse(y.pred == 1, 'Periodic', 'Linear')
  #table(pred, y)
  err <- mean(pred != y)
  return(err)
}


mserror <- c()
for(i in 1:500) {
  mserror <- rbind(mserror, missclassError(50,50))
}


summary(mserror)
boxplot(mserror, horizontal = T, xlab = 'Missclassification error rate')
#save(mserror, file = 'mserror.Rdata')
