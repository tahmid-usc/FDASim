rm(list = ls())
library(MASS)
library(kernlab)
library(mvtnorm)
library(Matrix)
library(optimx)
library(fdapace)
library(dplyr)
library(smoothmest)

source('code/RBF.R')
source('code/functionsNew.R')
#source('code/laplace.R')
#source('code/matern52.R')
#source('code/matern32.R')
#source('code/functions.R')
source('code/multistart.R')


# Plot mean function

# Mean function 
muf <- function(x) {
  fx <- .25 * dnorm(x, mean = -.95, sd = .06) + 
        .25 * dnorm(x, mean = -.5, sd = .05) + 
        .25 * dnorm(x, mean = 0, sd = .03) +   
        .25 * dnorm(x, mean = .5,sd = .02) 
  return(fx/max(fx))
}


muf <- function(x) {
  fx <- .25 * ddoublex(x, -.95,  .06) + 
    .25 * ddoublex(x,  -.5, .05) + 
    .25 * ddoublex(x, 0,  .03) +   
    .25 * ddoublex(x, .5, .02) 
  return(fx/max(fx))
}


curve(muf, -1, 1, lwd = 2)

#--------------------------------------------------

funcgen <- function(muf, theta) {
  x <- seq(-1,1,length.out = 100)
  mu <- muf(x)
  y <- mu + mvrnorm(1, rep(0,100), ker(x, l = theta[1], sigf = theta[2])) 
  y <- (y - mean(mu)) 
  return(list(x=x, y = y))
}


func <- funcgen(muf, theta = c(.01,.1))
plot(func$x, func$y, type = 'l', lwd = 3, 
     main = 'Simulated random functions', xlab = 'Time', ylab = 'Y',
     cex.main = 1.5, cex.lab = 1.5, ylim=c(-2, 2))
for(i in 1:10) {
  func <- funcgen(muf, theta = c(.1,.1))
  lines(func$x, func$y, lwd = 3)
}


#-------------------------------------------


fdata <- fdagen(n = 50, maxt = 10, muf = muf, theta = c(.01,.1,.01))

#plotting FDA
t <- seq(-1,1, length.out = 100)
#plot(t, muf(t), type = 'l', lwd = 4, ylim = c(-3,3))
plot(1, type="n", xlim=c(-1, 1), ylim=c(-2, 2), 
     main = 'Simulated functional data', xlab = 'Time', ylab = 'Y',
     cex.main = 1.5, cex.lab = 1.5)
for(i in 1:20){
  lines(fdata[fdata$id == i,1], fdata[fdata$id == i,2], type = 'b', lwd = 3)
}


fet <- feature(fdata)

time <- seq(-1, 1, length.out = 100)

posMean <- gpsmooth(time, fet)
plot(time, posMean, type = 'l', lwd = 3, ylim=c(-2, 2))
points(fdata$x, fdata$y, pch = 16)

#fpca

fpca.dt <- MakeFPCAInputs(IDs = fdata$id, fdata$x, fdata$y, sort = FALSE)
fpca <- FPCA(fpca.dt$Ly, fpca.dt$Lt, optns = list(dataType = 'Sparse'))

fdata <- arrange(fdata, x)
mu.fpca <- Lwls1D(fpca$bwMu, xin = fdata$x, yin = fdata$y, xout = time, kernel_type = 'gauss')

plot(time, mu.fpca, type = 'l', lwd = 3, ylim = c(-2,2))
points(fdata$x, fdata$y, pch = 16)

fmu <- fpcamu(time, fpca)


integrate(residfunc, lower = 0 , upper = 1, muf = muf, est = fpcamu, est_arg = fpca, subdivisions = 1000)
integrate(residfunc, lower = 0 , upper = 1, muf = muf, est = gpsmooth, est_arg = fet, subdivisions = 1000)


