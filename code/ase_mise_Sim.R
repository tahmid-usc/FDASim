# Replicate simulation and comparision


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


ASE <- function() {
  
  time <- seq(-1, 1, length.out = 100)
  # Generate training data----
  
  fdata <- fdagen(n = 20, maxt = 10, muf = muf1, theta = c(.01,.1,.01))
  # Fit FPCA
  
  fpca.dt <- MakeFPCAInputs(IDs = fdata$id, fdata$x, fdata$y, sort = FALSE)
  fpca <- FPCA(fpca.dt$Ly, fpca.dt$Lt, optns = list(dataType = 'Sparse'))
  
  #GP
  
  #source('code/laplace.R')
  #source('code/matern52.R')
  #source('code/matern32.R')
  
  fet <- feature(fdata)
  
  # Compute MISE
  
  mise.gp <- hcubature(residfunc, lower = -1 , upper = 1, muf = muf1, est = gpsmooth, est_arg = fet)$integral
  mise.fpca <- hcubature(residfunc, lower = -1 , upper = 1, muf = muf1, est = fpcamu, est_arg = fpca)$integral
  ase.fpca <- mean((residfunc(time, muf = muf1, est = fpcamu, est_arg = fpca)))
  ase.gp <- mean((residfunc(time, muf = muf1, est = gpsmooth, est_arg = fet)))
  
  return(cbind(mise.fpca, mise.gp, ase.fpca, ase.gp))
}

ase <- c()
for (i in 1:500) {
  ase <- rbind(ase, ASE())
}


apply(ase, 2, mean)

boxplot(ase[,1], ase[,2], horizontal = T, names = c('FPCA', 'GP'), cex.main = 1, 
        cex.lab = 1.5,
        main = 'Mean integrated square error (MISE) for estimating linear mean function')

#save(ase, file = 'ase_lin_m32.Rdata')

