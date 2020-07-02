rm(list = ls())
library(MASS)
library(kernlab)
library(mvtnorm)
library(Matrix)
library(optimx)
library(fdapace)
library(dplyr)
library(cubature)
library(smoothmest)

source('code/RBF.R')
source('code/functionsNew.R')
#source('code/laplace.R')
#source('code/matern52.R')
#source('code/matern32.R')
#source('code/functions.R')
source('code/multistart.R')


# Visualize mean function

t <- seq(0,1, length.out = 100)
curve(t,mu, 0, 1)

# Plot mean function

curve(muf1, -1, 1, lwd = 2)

#--------------------------------------------------

func <- funcgen(muf2, theta = c(.01,.1))
plot(func$x, func$y, type = 'l', lwd = 3, 
     main = 'Simulated random functions', xlab = 'Time', ylab = 'Y',
     cex.main = 1.5, cex.lab = 1.5, ylim=c(-2, 2))
for(i in 1:10) {
  func <- funcgen(muf1, theta = c(.1,.1))
  lines(func$x, func$y, lwd = 3)
}


#-------------------------------------------


fdata <- fdagen(n = 20, maxt = 10, muf = muf1, theta = c(.01,.1,.01))

#plotting FDA
t <- seq(-1,1, length.out = 100)
#plot(t, muf1(t), type = 'l', lwd = 4, ylim = c(-3,3))
plot(1, type="n", xlim=c(-1, 1), ylim=c(-2, 2), 
     main = 'Simulated functional data', xlab = 'Time', ylab = 'Y',
     cex.main = 1.5, cex.lab = 1.5)
for(i in 1:20){
  lines(fdata[fdata$id == i,1], fdata[fdata$id == i,2], type = 'b', lwd = 3)
}

#source('code/laplace.R')
#source('code/matern52.R')
#source('code/matern32.R')

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


plot(time, mu.fpca, type = 'l', lwd = 3, ylim = c(-1,2), xlab = 'Time', ylab = 'Y')
points(fdata$x, fdata$y, pch = 16)
points(time, posMean, type = 'l', lwd = 3, col = 2)
lines(time, muf1(time), lwd = 3, col = 3)
legend("bottomright", c("PACE", "GP", "True"), col = 1:3, lty = 1, bty = 'n')


# Testing integration method

f <- function(x) abs(sin(x))
f <- Vectorize(f)
integrate(f, -1, 1)
cubintegrate(f, -1 , 1)
hcubature(f, -1, 1)
mean(f(runif(1000, -1 , 1)))*2
mean(f(seq(-1,1, length.out = 1000)))

# Compute error
integrate(residfunc, lower = -1 , upper = 1, muf = muf1, est = fpcamu, 
          est_arg = fpca, subdivisions = 1000, abs.tol = .02)
integrate(residfunc, lower = -1 , upper = 1, muf = muf1, est = gpsmooth, 
          est_arg = fet, subdivisions = 1000)


cubintegrate(residfunc, lower = -1 , upper = 1, method = "pcubature", muf = muf1, est = gpsmooth, est_arg = fet)
cubintegrate(residfunc, lower = -1 , upper = 1, method = "pcubature", muf = muf1, est = fpcamu, est_arg = fpca)

# investiage effect of grid size
time <- seq(-1, 1, length.out = 100)
mean((residfunc(time, muf = muf1, est = fpcamu, est_arg = fpca)))
mean((residfunc(time, muf = muf1, est = gpsmooth, est_arg = fet)))


rmse <- function(n, method) {
  if(method == 'mc') time <- runif(n, -1, 1) else time <- seq(-1, 1, length.out = n)
  pca <- mean((residfunc(time, muf = muf1, est = fpcamu, est_arg = fpca)))
  gp <- mean((residfunc(time, muf = muf1, est = gpsmooth, est_arg = fet)))
  return(cbind(pca, gp))
  
}

vrmse <- Vectorize(rmse, vectorize.args = 'n')
#vrmse(n = 1:10, method = 'uni')
#vrmse(n = 1:10, method = 'mc')

n <- c(10, 50, 100, 200, 500, 750, 1000, 1500, 2000, 5000, 10000, 50000)
mse.unif <- vrmse(n, method = 'unif')
mse.mc <- vrmse(n, method = 'mc') * 2


mse.unif
mse.mc

plot(n, mse.unif[1,], col = 1, pch = 16, ylim = c(min(mse.unif), max(mse.unif)), type = 'l', 
     xlab = 'Grid size', ylab = 'MSE', lwd = 3, main = 'Fixed grid')
lines(n, mse.unif[2,], col = 2, pch = 16, lwd = 3)
legend('bottomright', c('PACE', 'GP'), lty = 1, col = 1:2, bty = 'n')

plot(n,mse.mc[1,], col = 1, pch = 16, ylim = c(min(mse.mc), max(mse.mc)), type = 'l', 
     xlab = 'Grid size', ylab = 'MSE', lwd = 3, main = 'Random grid')
lines(n, mse.mc[2,], col = 2, pch = 16, lwd = 3)
legend('bottomright', c('PACE', 'GP'), lty = 1, col = 1:2, bty = 'n')
