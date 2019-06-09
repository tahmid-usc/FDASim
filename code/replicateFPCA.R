# Replicate simulation and comparision


rm(list = ls())
library(fdapace)
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


ASE <- function() {
  # Generate training data----
  source('code/RBF.R')
  train <- gen(n = 50, theta1 = c(1000,1,.1), theta2 = c(1000,1,.1), maxt = 5)
  train <- split(train, train$z)
  source('code/matern32.R')
  # standardize
  
  my <- mean(c(train$train1$y, train$train2$y))
  stdy <- sd(c(train$train1$y, train$train2$y))
  
  train$train1$y <- (train$train1$y - my) / stdy
  train$train2$y <- (train$train2$y - my) / stdy
  
  # Fit FPCA
  
  dt1 <- MakeFPCAInputs(IDs = train$train1$id, train$train1$x, train$train1$y, sort = FALSE)
  dt2 <- MakeFPCAInputs(IDs = train$train2$id, train$train2$x, train$train2$y, sort = FALSE)
  
  
  pca1 <- FPCA(dt1$Ly, dt1$Lt, optns = list(dataType = 'Sparse'))
  pca2 <- FPCA(dt2$Ly, dt2$Lt, optns = list(dataType = 'Sparse'))
  
  # Fit GP
  
  fet <- lapply(train, feature)
  
  # Compare
  
  a1 <- 2
  b1 <- .1
  a2 <- 2.5
  b2 <- 1
  c2 <- 0.06
  
  mu1 <- a1 + b1 * pca1$workGrid
  mu1 <- (mu1 - my) / stdy
  mu2 <- a2 + b2 * sin(pca2$workGrid/5) + c2 * pca2$workGrid 
  mu2 <- (mu2 - my) / stdy
  
  mu1.gp <- gpsmooth(pca1$workGrid, fet$train1)
  mu2.gp <- gpsmooth(pca2$workGrid, fet$train2)
  
  #FPCA
  fpca.lin <- mean((mu1 - pca1$mu)^2)
  fpca.per <- mean((mu2 - pca2$mu)^2)
  # GP
  gp.lin <- mean((mu1 - mu1.gp)^2)
  gp.per <- mean((mu2 - mu2.gp)^2)
  
  return(cbind(fpca.lin, fpca.per, gp.lin, gp.per))
}

ase <- c()
for (i in 1:500) {
  ase <- rbind(ase, ASE())
}

apply(ase, 2, mean)

save(ase, file = 'asem32.Rdata')
