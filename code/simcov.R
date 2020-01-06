rm(list = ls())
library(MASS)
library(kernlab)
library(mvtnorm)
library(Matrix)
library(optimx)

source('code/rbfCov.R')
source('code/simcovFunc.R')
#source('code/laplace.R')
#source('code/matern52.R')
#source('code/matern32.R')
#source('code/multistart.R')



func <- funcgen(muf, z = 0, theta = c(1,.1))
plot(func$x, func$gt, type = 'l', lwd = 3, ylim = c(-3,3), 
     main = 'Simulated random functions with covariate effect', xlab = 'Time', ylab = 'Y',
     cex.main = 1.5, cex.lab = 1.5)
for(i in 1:10) {
  if(i <= 5){
  func <- funcgen(muf, z = 0, theta = c(1,.1))
  lines(func$x, func$gt, lwd = 3)}
  func <- funcgen(muf, z = 1, theta = c(1,.1))
  lines(func$x, func$gt, lwd = 3, col = 2)
}

legend('bottomright', c('z = 0', 'z = 1'), col = 1:2, bty = 'n', lwd = 2, cex = 1.5)



fdata <- fdagen(n = 50, maxt = 10, muf = muf, theta = c(.1,.1,.01))

#plotting FDA
t <- seq(0,1, length.out = 100)
#plot(t, muf(t), type = 'l', lwd = 4, ylim = c(-3,3))
plot(1, type="n", xlim=c(0, 1), ylim=c(-3, 3), 
     main = 'Simulated functional data with binary covariate z', xlab = 'Time', ylab = 'Y',
     cex.main = 1.5, cex.lab = 1.5)
for(i in 1:20){
  lines(fdata[fdata$id == i,1], fdata[fdata$id == i,2], type = 'b', lwd = 3, col = fdata[fdata$id == i,3] + 1)
}

legend('bottomright', c('z = 0', 'z = 1'), col = 1:2, lty = 1, bty = 'n', cex = 1.6)



fet <- feature(fdata)

mu0 <- gpsmooth(x = t, z = 0, fet)
mu1 <- gpsmooth(x = t, z = 1, fet)


plot(fdata[fdata$z ==0,]$x, fdata[fdata$z ==0,]$y, col = 2)
points(fdata[fdata$z ==1,]$x, fdata[fdata$z ==1,]$y, col = 3)
lines(t, mu0)
lines(t, mu1)