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
  
  #FPCA
  fpca.lin <- mean((mu1 - pca1$mu)^2)
  fpca.per <- mean((mu2 - pca2$mu)^2)
  
  
#  mu1.gp <- gpsmooth(pca1$workGrid, fet$train1)
#  mu2.gp <- gpsmooth(pca2$workGrid, fet$train2)
  
  #GP
  
  source('code/rbfmatern52.R')
  fet.rbf <- lapply(train, feature)
  rbf.lin <- mean((mu1 - gpsmooth(pca1$workGrid, fet.rbf$train1))^2)
  rbf.per <- mean((mu2 - gpsmooth(pca2$workGrid, fet.rbf$train2))^2)
  
  source('code/lapmatern52.R')
  fet.lap <- lapply(train, feature)
  lap.lin <- mean((mu1 - gpsmooth(pca1$workGrid, fet.lap$train1))^2)
  lap.per <- mean((mu2 - gpsmooth(pca2$workGrid, fet.lap$train2))^2)
  
  source('code/matern52.R')
  fet.m52 <- lapply(train, feature)
  m52.lin <- mean((mu1 - gpsmooth(pca1$workGrid, fet.m52$train1))^2)
  m52.per <- mean((mu2 - gpsmooth(pca2$workGrid, fet.m52$train2))^2)
  
  source('code/matern32matern52.R')
  fet.m32 <- lapply(train, feature)
  m32.lin <- mean((mu1 - gpsmooth(pca1$workGrid, fet.m32$train1))^2)
  m32.per <- mean((mu2 - gpsmooth(pca2$workGrid, fet.m32$train2))^2)
  

  # GP
  #gp.lin <- mean((mu1 - mu1.gp)^2)
  #gp.per <- mean((mu2 - mu2.gp)^2)
  
  return(cbind(fpca.lin, fpca.per, rbf.lin, lap.lin, m52.lin, m32.lin, 
               rbf.per, lap.per, m52.per, m32.per))
}

ase <- c()
for (i in 1:500) {
  ase <- rbind(ase, ASE())
}

apply(ase, 2, mean)

save(ase, file = 'ase.Rdata')

boxplot(ase[,c(1,3:6)], pch = 16, cex = .5, horizontal = T, 
        names = c('FPCA', 'RBF', 'Laplace', 'Matern 5/2', 'Matern 3/2'), xlab = 'Mean ASE', 
        main = 'Linear')


boxplot(ase[,-c(1,3:6)], pch = 16, cex = .5, horizontal = T, 
        names = c('FPCA', 'RBF', 'Laplace', 'Matern 5/2', 'Matern 3/2'), xlab = 'Mean ASE', 
        main = 'Periodic')