rm(list = ls())
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
source('code/multistart.R')

# Generate training data----

train <- gen(n = 50, theta1 = c(1000,1,.1), theta2 = c(1000,1,.1), maxt = 5)

train <- split(train, train$z)

# Test data

test <- gen(n = 50, theta1 = c(1000,1,.1), theta2 = c(1000,.1,.1), maxt = 5)

# Visualize ------

curve(2 + .1 * x, 0 , 100, lwd=4, main = 'Linear Population')
for (i in 1:10) {
  lines(train$train1[train$train1$id == i,]$x, train$train1[train$train1$id == i,]$y, lwd = 2, col = 'grey', type = 'b', pch = 16)
}


curve(2.5 + 1 * sin(x/5) + .06 * x, 0 , 100, lwd=4, main = 'Periodic Population')
for (i in 51:60) {
  lines(train$train2[train$train2$id == i,]$x, train$train2[train$train2$id == i,]$y, lwd = 2, col = 'grey', type = 'b', pch = 16)
}


curve(2 + .1 * x, 0 , 100, lwd = 4)
lines(0:100, 2.5 + 1 * sin(0:100/5) + .06 * 0:100, type = 'l', lwd = 4, col = 2)
points(train$train1$x, train$train1$y, pch = 16, col = 1, cex = .5)
points(train$train2$x, train$train2$y, pch = 16, col = 2, cex = .5)

# standardize

mx <- mean(c(train$train1$x, train$train2$x))
stdx <- sd(c(train$train1$x, train$train2$x))


my <- mean(c(train$train1$y, train$train2$y))
stdy <- sd(c(train$train1$y, train$train2$y))


train$train1$x <- (train$train1$x - mx) / stdx
train$train1$y <- (train$train1$y - my) / stdy


train$train2$x <- (train$train2$x - mx) / stdx
train$train2$y <- (train$train2$y - my) / stdy


test$x <- (test$x - mx) / stdx
test$y <- (test$y - my) / stdy

#---------------------------------

fet <- lapply(train, feature)

# Plot fitted
ageseq <- seq(-3,3, length = 100)
ageseq <- seq(min(train$train1$x, train$train2$x), max(train$train1$x, train$train2$x), length.out = 100)

plot(train$train1$x, train$train1$y, pch = 16, col = 1, 
     cex=.5, xlab='x', ylab='y')
lines(ageseq, 2 + .1 * ageseq, type = 'l', lwd = 2, col = 1, lty = 2)
lines(ageseq, 2.5 + 1 * sin(ageseq/5) + .06 * ageseq, type = 'l', lwd = 2, col = 2, , lty = 2)
points(train$train2$x, train$train2$y, pch = 16, col = 2, cex=.5)
lines(ageseq, gpsmooth(ageseq, fet$train1), type='l', lwd = 2, col = 1)
lines(ageseq, gpsmooth(ageseq, fet$train2), type='l', lwd = 2, col = 2)


#------------------------------


log.prob <- c()

for(i in unique(test$id)) {
  
  fit1 <- fit.gp(fet$train1, testx = test[test$id == i, 1], testy = test[test$id == i, 2])
  fit2 <- fit.gp(fet$train2, testx = test[test$id == i, 1], testy = test[test$id == i, 2])
  log.prob <- rbind(log.prob, c(fit1, fit2))
  
}

y <- c()
for (i in unique(test$id)) {
  y <- rbind(y, as.character((test[test$id == i, 4])[1]))
}

y.pred <- log.prob[,1] > log.prob[,2]
pred <- ifelse(y.pred == 1, 'train1', 'train2')
table(pred, y)
err <- mean(pred != y)
err