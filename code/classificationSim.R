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


train <- createfdt(50)
test <- createfdt(50)


#--------------------------------

train <- split(train,train$z)
fet <- lapply(train, feature)

curve(sin(x * 10) + x, 0 , 1, lwd = 5, xlab = 'Time', ylab = 'Y', cex.lab = 1.5, cex.axis = 1.5, ylim = c(-1.5, 2.5))
curve( 0 + 1.2 * x, 0 , 1, add = T, lwd = 5, col = 2)
for(i in unique(train$Linear$id)) {
  lines(train$Linear[train$Linear$id == i, 1], train$Linear[train$Linear$id == i, 2], col = 2, type = 'b')
}
for(i in unique(train$Periodic$id)) {
  lines(train$Periodic[train$Periodic$id == i, 1], train$Periodic[train$Periodic$id == i, 2], col = 1, type = 'b')
}
legend('bottomright', c( 'Periodic', 'Linear'), col = 1:2, lty = 1, bty = 'n', lwd = 2)


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
table(pred, y)
err <- mean(pred != y)
err


plotpred <- c(which(pred == y & pred == 'Periodic')[1], which(pred == y & pred == 'Linear')[1],
              which(pred != y & pred == 'Periodic')[1], which(pred != y & pred == 'Linear')[1])
plotpred
# plot


curve(sin(x * 10) + x, 0 , 1, lwd = 5, xlab = 'Time', ylab = 'Y', cex.lab = 1.5, cex.axis = 1.5, ylim = c(-1, 2))
curve( 0 + 1.2 * x, 0 , 1, add = T, lwd = 5, col = 2)
for(i in 1:4) {
  lines(test[test$id == plotpred[i], 1], test[test$id == plotpred[i], 2], type = 'b', col = i + 2, lwd = 5) 
}
legend(.8,-.5, c('Correct periodic', 'Correct linear', 'Incorrect periodic', 'Incorrect linear'), col = 3:6, lty = 1, bty = 'n', lwd = 2)


#save.image('classSim')
       