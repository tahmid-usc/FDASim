rm(list = ls())
library(MASS)
library(kernlab)
library(mvtnorm)
library(Matrix)
library(optimx)

source('code/RBF.R')
#source('code/laplace.R')
#source('code/matern52.R')
#source('code/matern32.R')
#source('code/functions.R')
#source('code/multistart.R')

# Mean function depends on time and covariate z
muf <- function(x, z) {
  #  (z+.8) * sin(x*10 + z) + z/2
  sin(x*10 + z) + z
  
}

#--------------------------------------------------
funcgen <- function(muf, z, theta) {
  x <- seq(0,1,length.out = 100)
  mu <- muf(x, z)
  gt <- mu + mvrnorm(1, rep(0,100), ker(x, l = theta[1], sigf = theta[2])) 
  return(list(x=x, gt = gt))
}


func <- funcgen(muf, z = 0, theta = c(1,.1))

colF <- topo.colors(10, alpha = .5)
colF <- rep(colF, each = 5)

plot(func$x, func$gt, type = 'n', lwd = 3, ylim = c(-3,7), col = colF[1], 
     main = 'Simulated random functions with covariate effect', xlab = 'Time', ylab = 'Y',
     cex.main = 1.5, cex.lab = 1.5 )

zval <- rep(seq(0,3,length.out = 10), each = 5)
for(i in 1:50) {
  
    func <- funcgen(muf, z = zval[i], theta = c(1,.1))
    lines(func$x, func$gt, lwd = 3, col = colF[i])
}


legend('bottomright', c(as.character(unique(round(zval,2)))), col = unique(colF), 
       horiz = T, bty = 'n', lwd = 3 )

# Second version with differenet color scheme


func <- funcgen(muf, z = 0, theta = c(1,.1))

plot(0, 0, type = 'n', lwd = 3, ylim = c(-2,4), xlim = c(0,1),
     main = 'Simulated random functions with covariate effect', xlab = 'Time', ylab = 'Y',
     cex.main = 1.5, cex.lab = 1.5 )

zval <- rep(seq(0,3,length.out = 10), each = 5)
for(i in 1:30) {
  
  func <- funcgen(muf, z = zval[i], theta = c(1,.1))
  lines(func$x, func$gt, lwd = 3, col = rgb(i/30,0,0))
}

legend('bottomright', as.character(round(c(0, 1/3, 1), 3)), col = rgb(c(0, 1/3, 1),0,0), 
       bty = 'n', lwd = 5, cex = 1.5, horiz = T)


#-------------------


#-------------------------------------------

fdagen <- function(n = 100, maxt = 20, muf, theta = rep(1,3)) {
  
  train <- data.frame()
  #n number of functions in the sample data
  n.time <- sample(1:maxt, size=n, replace=T) 
  
  for(i in 1:n) {
    id <- rep(i, n.time[i])
    z <- rep(runif(1, 0, 3), n.time[i])
    x <- sort(runif(n.time[i], 0 , 1))
    mu <- muf(x,z)
    gt <- mvrnorm(1, rep(0, n.time[i]), ker(x, l = theta[1], sigf = theta[2])) 
    #y <- mu + gt
    y <- mu + gt + rnorm(n.time[i], 0, theta[3])
    train <- rbind(train, cbind(x, y, z, id))
  }
  return(train)
}



fdata <- fdagen(n = 50, maxt = 10, muf = muf, theta = c(.1,.1,.01))

#plotting FDA
t <- seq(0,1, length.out = 100)
#plot(t, muf(t), type = 'l', lwd = 4, ylim = c(-3,3))
plot(1, type="n", xlim=c(0, 1), ylim=c(-1, 4), 
     main = 'Simulated functional data with continuous covariate', xlab = 'Time', ylab = 'Y',
     cex.main = 1.5, cex.lab = 1.5)
for(i in 1:20){
  lines(fdata[fdata$id == i,1], fdata[fdata$id == i,2], type = 'b', lwd = 3, col = rgb(fdata[fdata$id == i,3]/3,0,0))
}

legend('bottomright', as.character(round(c(0, 1/3, 1), 3)), col = rgb(c(0, 1/3, 1),0,0), 
       bty = 'n', lwd = 5, cex = 1.5, horiz = T)

#---------------------

fet <- feature(fdata)

zval <- seq(0,3,length.out = 20)
mu <- c()
for(i in 1:20){
  mu <- cbind(mu, gpsmooth(x = t, z = zval[i], fet)) 
}


plot(fdata$x, fdata$y, pch = 16, col = rgb(fdata[,3]/3,0,0),
     main = 'Fitting GP model with continuous covariate', xlab = 'Time', ylab = 'Y', cex.lab = 1.5, cex.axis = 1.5)
for(i in 1:20){
  lines(t, mu[,i], type = 'l', lwd = 3, col = rgb(i/20,0,0))
}

legend('bottomright', as.character(round(c(0, 1/3, 1), 3)), col = rgb(c(0, 1/3, 1),0,0), 
       bty = 'n', lwd = 5, cex = 1.5, horiz = T)





