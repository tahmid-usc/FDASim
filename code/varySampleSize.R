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


sampleSize = 50

# generate data

source('code/RBF.R')
fdata1 <- fdagen(n = sampleSize, maxt = 10, muf = muf1, theta = c(1, 5, 2))

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
ase.muf1 <- ase(N, method = 'unif', muf = muf1, est = gpsmooth, est_arg = fet.muf1.rbf)

source('code/laplace.R')
ase.muf1 <- cbind(ase.muf1, ase(N, method = 'unif', muf = muf1, est = gpsmooth, est_arg = fet.muf1.lap))

source('code/matern52.R')
ase.muf1 <- cbind(ase.muf1, ase(N, method = 'unif', muf = muf1, est = gpsmooth, est_arg = fet.muf1.m52))

source('code/matern32.R')
ase.muf1 <- cbind(ase.muf1, ase(N, method = 'unif', muf = muf1, est = gpsmooth, est_arg = fet.muf1.m32))

ase.muf1 <- cbind(ase.muf1, ase(N, method = 'unif', muf = muf1, est = fpcamu, est_arg = fpca1))

sase.muf1 <- ase.muf1/var(fdata1$y)


kername <- c('Gaussian', 'Laplace', 'Matern_5/2', 'Matern 3/2', 'PACE')
SASE <- round(rbind(sase.muf1), 4)  
SASE <- cbind(c('Mix of normal'), SASE)
colnames(SASE) <- c('Method', kername)


# Mean Integrated Squared Error (MISE)

# # Base R
# source('code/RBF.R')
# rint.muf1 <- possibly(integrate(residfunc, lower = 0 , upper = 1, muf = muf1, est = gpsmooth, 
#                        est_arg = fet.muf1.rbf, subdivisions = 10000)$value, otherwise = NA)
# source('code/laplace.R')
# rint.muf1 <- cbind(rint.muf1,
#                    possibly(integrate(residfunc, lower = 0 , upper = 1, muf = muf1, est = gpsmooth, 
#                              est_arg = fet.muf1.lap, subdivisions = 10000)$value, otherwise = NA))
# source('code/matern52.R')
# rint.muf1 <- cbind(rint.muf1,
#                    possibly(integrate(residfunc, lower = 0 , upper = 1, muf = muf1, est = gpsmooth, 
#                              est_arg = fet.muf1.m52, subdivisions = 10000)$value, otherwise = NA))
# source('code/matern32.R')
# rint.muf1 <- cbind(rint.muf1,
#                    possibly(integrate(residfunc, lower = 0 , upper = 1, muf = muf1, est = gpsmooth, 
#                              est_arg = fet.muf1.m32, subdivisions = 10000)$value, otherwise = NA))
# rint.muf1 <- cbind(rint.muf1,
#                    possibly(integrate(residfunc, lower = 0, upper = 1, muf = muf1, est = fpcamu, 
#                              est_arg = fpca1, subdivisions = 10000)$value, otherwise = NA))
# # standardizing
# rint.muf1 <- rint.muf1/var(fdata1$y)
# 
# 
# RINT <- round(rbind(rint.muf1), 4)  
# RINT <- cbind(c('Mix of normal'), RINT)
# colnames(RINT) <- c('Method', kername)
# 

# Cubature

source('code/RBF.R')
cubint.muf1 <- cubintegrate(residfunc, lower = 0 , upper = 1, method = "pcubature", muf = muf1, est = gpsmooth, est_arg = fet.muf1.rbf)$integral
source('code/laplace.R')
cubint.muf1 <- cbind(cubint.muf1,
                     cubintegrate(residfunc, lower = 0 , upper = 1, method = "pcubature", muf = muf1, est = gpsmooth, est_arg = fet.muf1.lap)$integral)
source('code/matern52.R')
cubint.muf1 <- cbind(cubint.muf1,
                     cubintegrate(residfunc, lower = 0 , upper = 1, method = "pcubature", muf = muf1, est = gpsmooth, est_arg = fet.muf1.m52)$integral)
source('code/matern32.R')
cubint.muf1 <- cbind(cubint.muf1,
                     cubintegrate(residfunc, lower = 0 , upper = 1, method = "pcubature", muf = muf1, est = gpsmooth, est_arg = fet.muf1.m32)$integral)
cubint.muf1 <- cbind(cubint.muf1,
                     cubintegrate(residfunc, lower = 0 , upper = 1, method = "pcubature", muf = muf1, est = fpcamu, est_arg = fpca1)$integral)

#standardize
cubint.muf1 <- cubint.muf1/var(fdata1$y)

CUBINT <- round(rbind(cubint.muf1), 4)  
CUBINT <- cbind(c('ISE: Mix of normal'), CUBINT)
colnames(CUBINT) <- c('Method', kername)



# MC method

source('code/RBF.R')
ise.muf1 <- ase(N, method = 'mc', muf = muf1, est = gpsmooth, est_arg = fet.muf1.rbf) 
source('code/laplace.R')
ise.muf1 <- cbind(ise.muf1, ase(N, method = 'mc', muf = muf1, est = gpsmooth, est_arg = fet.muf1.lap)) 
source('code/matern52.R')
ise.muf1 <- cbind(ise.muf1, ase(N, method = 'mc', muf = muf1, est = gpsmooth, est_arg = fet.muf1.m52)) 
source('code/matern32.R')
ise.muf1 <- cbind(ise.muf1, ase(N, method = 'mc', muf = muf1, est = gpsmooth, est_arg = fet.muf1.m32)) 

ise.muf1 <- cbind(ise.muf1, ase(N, method = 'mc', muf = muf1, est = fpcamu, est_arg = fpca1)) 

#standardizing
ise.muf1 <- ise.muf1 / var(fdata1$y)

MC <- round(rbind(ise.muf1), 4)  
MC <- cbind(c('ISE: Mix of normal'), MC)
colnames(MC) <- c('Method', kername)


Result <- round(rbind(sase.muf1, cubint.muf1, ise.muf1), 4)  
Result <- cbind(c('SASE', 'ISE: Cubature', 'ISE: MC'), Result)
colnames(Result) <- c('Method', kername)

list(SASE= SASE, CUBINT = CUBINT, MC = MC, Result = Result)




#Result10 <- Result
#Result20 <- Result
#Result30 <- Result
#Result40 <- Result
#Result50 <- Result

#save(list = ls(all.names = TRUE), file = "sampleSizeVary.RData", envir = .GlobalEnv)
#save(resultMat, file = 'varySamMat.Rdata')

par(mfrow = c(1,3))

n.vec <- c(10,20,30,40, 50)
resultMat <- rbind(Result10,Result20,Result30,Result40,Result50)
resultMat <- cbind(rep(n.vec, each = 3), resultMat)

par(mfrow = c(1,3))
plot(n.vec, c(Result10[1,2],Result20[1,2],Result30[1,2],Result40[1,2],Result50[1,2]), type = 'l', ylim = c(0,1))
lines(n.vec, c(Result10[1,3],Result20[1,3],Result30[1,3],Result40[1,3],Result50[1,3]), col = 2)
lines(n.vec, c(Result10[1,4],Result20[1,4],Result30[1,4],Result40[1,4],Result50[1,4]), col = 3)
lines(n.vec, c(Result10[1,5],Result20[1,5],Result30[1,5],Result40[1,5],Result50[1,5]), col = 4)
lines(n.vec, c(Result10[1,6],Result20[1,6],Result30[1,6],Result40[1,6],Result50[1,6]), col = 5)


plot(n.vec, c(Result10[2,2],Result20[2,2],Result30[2,2],Result40[2,2],Result50[2,2]), type = 'l', ylim = c(0,1))
lines(n.vec, c(Result10[2,3],Result20[2,3],Result30[2,3],Result40[2,3],Result50[2,3]), col = 2)
lines(n.vec, c(Result10[2,4],Result20[2,4],Result30[2,4],Result40[2,4],Result50[2,4]), col = 3)
lines(n.vec, c(Result10[2,5],Result20[2,5],Result30[2,5],Result40[2,5],Result50[2,5]), col = 4)
lines(n.vec, c(Result10[2,6],Result20[2,6],Result30[2,6],Result40[2,6],Result50[2,6]), col = 5)


plot(n.vec,  c(Result10[3,2],Result20[3,2],Result30[3,2],Result40[3,2],Result50[3,2]), type = 'l', ylim = c(0,1))
lines(n.vec, c(Result10[3,3],Result20[3,3],Result30[3,3],Result40[3,3],Result50[3,3]), col = 2)
lines(n.vec, c(Result10[3,4],Result20[3,4],Result30[3,4],Result40[3,4],Result50[3,4]), col = 3)
lines(n.vec, c(Result10[3,5],Result20[3,5],Result30[3,5],Result40[3,5],Result50[3,5]), col = 4)
lines(n.vec, c(Result10[3,6],Result20[3,6],Result30[3,6],Result40[3,6],Result50[3,6]), col = 5)


plot(n.vec, resultMat[seq(1, 13, by = 3),3], 
     type = 'b', xlab = 'Sample size', ylab = 'SASE', cex.lab = 1.5, lwd = 3, 
     ylim = c(0, .5))
for (i in 4:7) lines(n.vec, resultMat[seq(1, 13, by = 3),i], lwd = 3, col = i - 2)
legend('topright', c('Gaussian', 'Laplace', 'Matern 5/2', 'Matern 3/2', 'PACE'), 
       lty = 1, bty = 'n', col = 1:5, lwd = 3)

plot(n.vec, resultMat[seq(2, 14, by = 3),3], 
     type = 'b', xlab = 'Sample size', ylab = 'ISE: Cubature', cex.lab = 1.5, lwd = 3, 
     ylim = c(0, .5))
for (i in 4:7) lines(n.vec, resultMat[seq(2, 14, by = 3),i], lwd = 3, col = i - 2)
legend('topright', c('Gaussian', 'Laplace', 'Matern 5/2', 'Matern 3/2', 'PACE'), 
       lty = 1, bty = 'n', col = 1:5, lwd = 3)


plot(n.vec, resultMat[seq(3, 15, by = 3),3], 
     type = 'b', xlab = 'Sample size', ylab = 'ISE: MC', cex.lab = 1.5, lwd = 3, 
     ylim = c(0, 1.2))
for (i in 4:7) lines(n.vec, resultMat[seq(3, 15, by = 3),i], lwd = 3, col = i - 2)
legend('topright', c('Gaussian', 'Laplace', 'Matern 5/2', 'Matern 3/2', 'PACE'), 
       lty = 1, bty = 'n', col = 1:5, lwd = 3)

