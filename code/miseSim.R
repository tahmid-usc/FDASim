# Replicate simulation and comparision


rm(list = ls())
library(fdapace)
library(MASS)
library(kernlab)
library(mvtnorm)
library(Matrix)
library(optimx)
library(dplyr)

source('code/RBF.R')
#source('code/laplace.R')
#source('code/matern52.R')
#source('code/matern32.R')
source('code/functionsNew.R')


ASE <- function() {
  
  time <- seq(0, 1, length.out = 100)
  # Generate training data----
  
  fdata <- fdagen(n = 50, maxt = 10, muf = muf, theta = c(.1,.1,.01))
  
  # Fit FPCA
  
  fpca.dt <- MakeFPCAInputs(IDs = fdata$id, fdata$x, fdata$y, sort = FALSE)
  fpca <- FPCA(fpca.dt$Ly, fpca.dt$Lt, optns = list(dataType = 'Sparse'))

  #GP
  
  #source('code/laplace.R')
  #source('code/matern52.R')
  source('code/matern32.R')
  
  fet <- feature(fdata)
  
  # Compute MISE
  
  mise.fpca <- integrate(residfunc, lower = 0 , upper = 1, muf = muf, est = fpcamu, est_arg = fpca, subdivisions = 10000)$value
  mise.gp <-integrate(residfunc, lower = 0 , upper = 1, muf = muf, est = gpsmooth, est_arg = fet, subdivisions = 10000)$value
  
  
  return(cbind(mise.fpca, mise.gp))
}

ase <- c()
for (i in 1:500) {
  ase <- rbind(ase, ASE())
}

apply(ase, 2, mean)

boxplot(ase[,1], ase[,2], horizontal = T, names = c('FPCA', 'GP'), cex.main = 1, 
        cex.lab = 1.5,
        main = 'Mean integrated square error (MISE) for estimating linear mean function')

save(ase, file = 'ase_lin_m32.Rdata')

