#loading packages
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
library(kableExtra)
library(colorspace)
library(purrr)

source('code/RBF.R')
source('code/functionsTest.R')
#source('code/functionsNew.R')
#source('code/laplace.R')
#source('code/matern52.R')
#source('code/matern32.R')
source('code/multistart.R')


errorEstMean <- function(sampleSize = 5, muf = muf1, maxt = 10) {

# generate data

source('code/RBF.R')
fdata1 <- fdagen(n = sampleSize, maxt = maxt, muf = muf, theta = c(1, 5, 2))

# GP fit for muf1

time <- seq(0, 1, length.out = 100)

#Gaussian kernel
source('code/RBF.R')
fet.muf1.rbf <- feature(fdata1)

#Laplace kernel
source('code/laplace.R')
fet.muf1.lap <- feature(fdata1)

# Mtern 5/2 kernel

source('code/matern52.R')
fet.muf1.m52 <- feature(fdata1)

# Matern 3/2
source('code/matern32.R')
fet.muf1.m32 <- feature(fdata1)


# PACE
fpca.dt1 <- MakeFPCAInputs(IDs = fdata1$id, fdata1$x, fdata1$y, sort = FALSE)
fpca1 <- FPCA(fpca.dt1$Ly, fpca.dt1$Lt, optns = list(dataType = 'Sparse'))


# Standardized Average Squared Error (SASE)

N <- 10000


source('code/RBF.R')
ase.muf1 <- ase(N, method = 'unif', muf = muf, est = gpsmooth, est_arg = fet.muf1.rbf)

source('code/laplace.R')
ase.muf1 <- cbind(ase.muf1, ase(N, method = 'unif', muf = muf, est = gpsmooth, est_arg = fet.muf1.lap))

source('code/matern52.R')
ase.muf1 <- cbind(ase.muf1, ase(N, method = 'unif', muf = muf, est = gpsmooth, est_arg = fet.muf1.m52))

source('code/matern32.R')
ase.muf1 <- cbind(ase.muf1, ase(N, method = 'unif', muf = muf, est = gpsmooth, est_arg = fet.muf1.m32))

ase.muf1 <- cbind(ase.muf1, ase(N, method = 'unif', muf = muf, est = fpcamu, est_arg = fpca1))

sase.muf1 <- ase.muf1/var(fdata1$y)


# Cubature

source('code/RBF.R')
cubint.muf1 <- cubintegrate(residfunc, lower = 0 , upper = 1, method = "pcubature", muf = muf, est = gpsmooth, est_arg = fet.muf1.rbf)$integral
source('code/laplace.R')
cubint.muf1 <- cbind(cubint.muf1,
                     cubintegrate(residfunc, lower = 0 , upper = 1, method = "pcubature", muf = muf, est = gpsmooth, est_arg = fet.muf1.lap)$integral)
source('code/matern52.R')
cubint.muf1 <- cbind(cubint.muf1,
                     cubintegrate(residfunc, lower = 0 , upper = 1, method = "pcubature", muf = muf, est = gpsmooth, est_arg = fet.muf1.m52)$integral)
source('code/matern32.R')
cubint.muf1 <- cbind(cubint.muf1,
                     cubintegrate(residfunc, lower = 0 , upper = 1, method = "pcubature", muf = muf, est = gpsmooth, est_arg = fet.muf1.m32)$integral)
cubint.muf1 <- cbind(cubint.muf1,
                     cubintegrate(residfunc, lower = 0 , upper = 1, method = "pcubature", muf = muf, est = fpcamu, est_arg = fpca1)$integral)


# MC method

source('code/RBF.R')
ise.muf1 <- ase(N, method = 'mc', muf = muf, est = gpsmooth, est_arg = fet.muf1.rbf) 
source('code/laplace.R')
ise.muf1 <- cbind(ise.muf1, ase(N, method = 'mc', muf = muf, est = gpsmooth, est_arg = fet.muf1.lap)) 
source('code/matern52.R')
ise.muf1 <- cbind(ise.muf1, ase(N, method = 'mc', muf = muf, est = gpsmooth, est_arg = fet.muf1.m52)) 
source('code/matern32.R')
ise.muf1 <- cbind(ise.muf1, ase(N, method = 'mc', muf = muf, est = gpsmooth, est_arg = fet.muf1.m32)) 

ise.muf1 <- cbind(ise.muf1, ase(N, method = 'mc', muf = muf, est = fpcamu, est_arg = fpca1)) 


return(list(SASE= sase.muf1, CUBINT = cubint.muf1, MC = ise.muf1))

}

sampleSize = 50
muf = muf1
maxt = 10

result <- c()
for(i in 1:100) {
  err <- errorEstMean(sampleSize, muf, maxt)
  result$SASE <- rbind(result$SASE, err$SASE) 
  result$CUBINT <- rbind(result$CUBINT, err$CUBINT)
  result$MC <- rbind(result$MC, err$MC)
}


#save(list = ls(all.names = TRUE), file = "result/mse_muf1_n50_maxt10.RData", envir = .GlobalEnv)
