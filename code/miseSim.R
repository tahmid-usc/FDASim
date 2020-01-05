# Replicate simulation and comparision


rm(list = ls())
library(fdapace)
library(MASS)
library(kernlab)
library(mvtnorm)
library(Matrix)
library(optimx)

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

save(ase, file = 'ase_per.Rdata')

boxplot(ase[,c(1,3:6)], pch = 16, cex = .5, horizontal = T, 
        names = c('FPCA', 'RBF', 'Laplace', 'Matern 5/2', 'Matern 3/2'), xlab = 'Mean ASE', 
        main = 'Linear')


boxplot(ase[,-c(1,3:6)], pch = 16, cex = .5, horizontal = T, 
        names = c('FPCA', 'RBF', 'Laplace', 'Matern 5/2', 'Matern 3/2'), xlab = 'Mean ASE', 
        main = 'Periodic')