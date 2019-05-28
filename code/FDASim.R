rm(list = ls())
library(MASS)
library(kernlab)
library(mvtnorm)
library(Matrix)
library(optimx)

source('code/RBF.R')
source('code/matern52.R')
source('code/functions.R')
source('code/multistart.R')

#Generate training data----

train <- gen(n=50, theta = c(1000,.5,.2), maxt = 10)

train <- split(train, train$z)

# Test data

test <- gen(n=50, theta = c(1000,.5,.2), maxt = 10)

# Visualize ------

curve(2 + .1 * x, 0 , 100, lwd=4, main = 'Linear Population')
for (i in 1:10) {
  lines(train1[train1$id == i,]$x, train1[train1$id == i,]$y, lwd = 2, col = 'grey', type = 'b', pch = 16)
}


curve(2.5 + 1 * sin(x/5) + .06 * x, 0 , 100, lwd=4, main = 'Periodic Population')
for (i in 1:10) {
  lines(train2[train2$id == i,]$x, train2[train2$id == i,]$y, lwd = 2, col = 'grey', type = 'b', pch = 16)
}


curve(2 + .1 * x, 0 , 100, lwd = 4)
lines(0:100, 2.5 + 1 * sin(0:100/5) + .06 * 0:100, type = 'l', lwd = 4, col = 2)
points(train1$x, train1$y, pch = 16, col = 1, cex = .5)
points(train2$x, train2$y, pch = 16, col = 2, cex = .5)

# standardize

mx <- mean(c(train1$x, train2$x))
stdx <- sd(c(train1$x, train2$x))


my <- mean(c(train1$y, train2$y))
stdy <- sd(c(train1$y, train2$y))


train1$x <- (train1$x - mx) / stdx
train1$y <- (train1$y - my) / stdy


train2$x <- (train2$x - mx) / stdx
train2$y <- (train2$y - my) / stdy


train <- list(train1 = train1, train2 = train2)

test$x <- (test$x - mx) / stdx
test$y <- (test$y - my) / stdy

#---------------------------------

fet <- lapply(train, feature)

# Plot fitted
ageseq <- seq(0,100, length = 100)

plot(fet$train1$trainx, fet$train1$trainy, pch = 16, col = 1, 
     cex=.5, xlab='Age', ylab='Spinal bone mineral density')
points(fet$train2$trainx, fet$train2$trainy, pch = 16, col = 2, cex=.5)
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