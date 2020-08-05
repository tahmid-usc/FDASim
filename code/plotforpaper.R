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

source('./code/RBF.R')
source('./code/functionsTest.R')
#source('../code/functionsNew.R')
#source('code/laplace.R')
#source('code/matern52.R')
#source('code/matern32.R')
source('./code/multistart.R')

#load previous results
load("E:/FDASim/markdown/.RData")




pdf(file = 'plot_paper/mu1.pdf', width = 14, height = 6)
plot(1, type="n", xlim=c(0, 1), ylim=c(-10, 28), xlab = 'Time', ylab = 'Y',
     cex.main = 1.5, cex.lab = 1.5)
for(i in 1:20){
  lines(fdata1[fdata1$id == i,1], fdata1[fdata1$id == i,2], type = 'b', lwd = 2, col = rgb(0,0,0,.5))
}
lines(t, muf1(t), type = 'l', lwd = 4)
dev.off()




pdf(file = 'plot_paper/mu2.pdf', width = 14, height = 6)
plot(1, type="n", xlim=c(0, 1), ylim=c(-10, 28), xlab = 'Time', ylab = 'Y',
     cex.main = 1.5, cex.lab = 1.5)
for(i in 1:20){
  lines(fdata2[fdata2$id == i,1], fdata2[fdata2$id == i,2], type = 'b', lwd = 2, col = rgb(0,0,0,.5))
}
lines(t, muf2(t), type = 'l', lwd = 4)
dev.off()






pdf(file = 'plot_paper/mu1_fit.pdf', width = 14, height = 6)
plot(time, muf1(time), type = 'l', lwd = 3, ylim = c(-5,28), xlab = 'Time', ylab = 'Y', cex.lab = 1.5, col = 1)
points(fdata1$x, fdata1$y, pch = 16, cex = .8)
lines(time, posMean.muf1.rbf, type = 'l', lwd = 3, col = 2)
lines(time, posMean.muf1.lap, type = 'l', lwd = 3, col = 3)
lines(time, posMean.muf1.m52, type = 'l', lwd = 3, col = 4)
lines(time, posMean.muf1.m32, type = 'l', lwd = 3, col = 5)
lines(time, fmu1, lwd = 3, col = 6)
legend("topright", c("True", "Gaussian", 'Laplace', 'Matern 5/2', 'Matern 3/2', "PACE"), col = 1:6, lty = 1, bty = 'n', horiz = F, lwd = 3)
dev.off()






pdf(file = 'plot_paper/mu2_fit.pdf', width = 14, height = 6)
plot(time, fmu2, type = 'l', lwd = 3, ylim = c(-5,28), xlab = 'Time', cex.lab = 1.5, ylab = 'Y', col = 6)
points(fdata2$x, fdata2$y, pch = 16, cex = .8)
lines(time, posMean.muf2.rbf, type = 'l', lwd = 3, col = 2)
lines(time, posMean.muf2.lap, type = 'l', lwd = 3, col = 3)
lines(time, posMean.muf2.m52, type = 'l', lwd = 3, col = 4)
lines(time, posMean.muf2.m32, type = 'l', lwd = 3, col = 5)
lines(time, muf2(time), lwd = 3, col = 1)
legend("topright", c("True", "Gaussian", 'Laplace', 'Matern 5/2', 'Matern 3/2', "PACE"), col = 1:6, lty = 1, bty = 'n', horiz = F, lwd = 3)
dev.off()














pdf(file = 'plot_paper/credible_band_mu1_rbf.pdf', width = 14, height = 6)
source('code/RBF.R')
post.dist <- credible_band(t, fet.muf1.rbf)
ll <- post.dist$pred - 1.96 * sqrt(diag(post.dist$sigma))
ul <- post.dist$pred + 1.96 * sqrt(diag(post.dist$sigma))

plot(time, muf1(time), type = 'l', lwd = 3, ylim = c(-5,28), xlab = 'Time', ylab = 'Y', cex.lab = 1.5, col = 1)
points(fdata1$x, fdata1$y, pch = 16, cex = .8)
lines(time, post.dist$pred, type = 'l', lwd = 3, col = 2)
lines(time, ll, type = 'l', lwd = 3, col = 2, lty = 2)
lines(time, ul, type = 'l', lwd = 3, col = 2, lty = 2)
legend("topright", c("True", "Gaussian"), col = 1:2, lty = 1, bty = 'n', horiz = F, lwd = 3)
dev.off()

pdf(file = 'plot_paper/credible_band_mu1_laplace.pdf', width = 14, height = 6)
source('code/laplace.R')
post.dist <- credible_band(t, fet.muf1.lap)
ll <- post.dist$pred - 1.96 * sqrt(diag(post.dist$sigma))
ul <- post.dist$pred + 1.96 * sqrt(diag(post.dist$sigma))

plot(time, muf1(time), type = 'l', lwd = 3, ylim = c(-5,28), xlab = 'Time', ylab = 'Y', cex.lab = 1.5, col = 1)
points(fdata1$x, fdata1$y, pch = 16, cex = .8)
lines(time, post.dist$pred, type = 'l', lwd = 3, col = 3)
lines(time, ll, type = 'l', lwd = 3, col = 3, lty = 2)
lines(time, ul, type = 'l', lwd = 3, col = 3, lty = 2)
legend("topright", c("True", "Laplace"), col = c(1,3), lty = 1, bty = 'n', horiz = F, lwd = 3)
dev.off()

pdf(file = 'plot_paper/credible_band_mu1_m52.pdf', width = 14, height = 6)
source('code/matern52.R')
post.dist <- credible_band(t, fet.muf1.m52)
ll <- post.dist$pred - 1.96 * sqrt(diag(post.dist$sigma))
ul <- post.dist$pred + 1.96 * sqrt(diag(post.dist$sigma))

plot(time, muf1(time), type = 'l', lwd = 3, ylim = c(-5,28), xlab = 'Time', ylab = 'Y', cex.lab = 1.5, col = 1)
points(fdata1$x, fdata1$y, pch = 16, cex = .8)
lines(time, post.dist$pred, type = 'l', lwd = 3, col = 4)
lines(time, ll, type = 'l', lwd = 3, col = 4, lty = 2)
lines(time, ul, type = 'l', lwd = 3, col = 4, lty = 2)
legend("topright", c("True", "Matern 5/2"), col = c(1,4), lty = 1, bty = 'n', horiz = F, lwd = 3)
dev.off()

pdf(file = 'plot_paper/credible_band_mu1_m32.pdf', width = 14, height = 6)
source('code/matern32.R')
post.dist <- credible_band(t, fet.muf1.m32)
ll <- post.dist$pred - 1.96 * sqrt(diag(post.dist$sigma))
ul <- post.dist$pred + 1.96 * sqrt(diag(post.dist$sigma))

plot(time, muf1(time), type = 'l', lwd = 3, ylim = c(-5,28), xlab = 'Time', ylab = 'Y', cex.lab = 1.5, col = 1)
points(fdata1$x, fdata1$y, pch = 16, cex = .8)
lines(time, post.dist$pred, type = 'l', lwd = 3, col = 5)
lines(time, ll, type = 'l', lwd = 3, col = 5, lty = 2)
lines(time, ul, type = 'l', lwd = 3, col = 5, lty = 2)
legend("topright", c("True", "Matern 3/2"), col = c(1,5), lty = 1, bty = 'n', horiz = F, lwd = 3)
dev.off()





pdf(file = 'plot_paper/credible_band_mu2_rbf.pdf', width = 14, height = 6)
source('code/RBF.R')
post.dist <- credible_band(t, fet.muf2.rbf)
ll <- post.dist$pred - 1.96 * sqrt(diag(post.dist$sigma))
ul <- post.dist$pred + 1.96 * sqrt(diag(post.dist$sigma))

plot(time, muf2(time), type = 'l', lwd = 3, ylim = c(-5,28), xlab = 'Time', ylab = 'Y', cex.lab = 1.5, col = 1)
points(fdata2$x, fdata2$y, pch = 16, cex = .8)
lines(time, post.dist$pred, type = 'l', lwd = 3, col = 2)
lines(time, ll, type = 'l', lwd = 3, col = 2, lty = 2)
lines(time, ul, type = 'l', lwd = 3, col = 2, lty = 2)
legend("topright", c("True", "Gaussian"), col = 1:2, lty = 1, bty = 'n', horiz = F, lwd = 3)
dev.off()

pdf(file = 'plot_paper/credible_band_mu2_laplace.pdf', width = 14, height = 6)
source('code/laplace.R')
post.dist <- credible_band(t, fet.muf2.lap)
ll <- post.dist$pred - 1.96 * sqrt(diag(post.dist$sigma))
ul <- post.dist$pred + 1.96 * sqrt(diag(post.dist$sigma))

plot(time, muf2(time), type = 'l', lwd = 3, ylim = c(-5,28), xlab = 'Time', ylab = 'Y', cex.lab = 1.5, col = 1)
points(fdata2$x, fdata2$y, pch = 16, cex = .8)
lines(time, post.dist$pred, type = 'l', lwd = 3, col = 3)
lines(time, ll, type = 'l', lwd = 3, col = 3, lty = 2)
lines(time, ul, type = 'l', lwd = 3, col = 3, lty = 2)
legend("topright", c("True", "Laplace"), col = c(1,3), lty = 1, bty = 'n', horiz = F, lwd = 3)
dev.off()



pdf(file = 'plot_paper/credible_band_mu2_m52.pdf', width = 14, height = 6)
source('code/matern52.R')
post.dist <- credible_band(t, fet.muf2.m52)
ll <- post.dist$pred - 1.96 * sqrt(diag(post.dist$sigma))
ul <- post.dist$pred + 1.96 * sqrt(diag(post.dist$sigma))

plot(time, muf2(time), type = 'l', lwd = 3, ylim = c(-5,28), xlab = 'Time', ylab = 'Y', cex.lab = 1.5, col = 1)
points(fdata2$x, fdata2$y, pch = 16, cex = .8)
lines(time, post.dist$pred, type = 'l', lwd = 3, col = 4)
lines(time, ll, type = 'l', lwd = 3, col = 4, lty = 2)
lines(time, ul, type = 'l', lwd = 3, col = 4, lty = 2)
legend("topright", c("True", "Matern 5/2"), col = c(1,4), lty = 1, bty = 'n', horiz = F, lwd = 3)
dev.off()


pdf(file = 'plot_paper/credible_band_mu2_m32.pdf', width = 14, height = 6)
source('code/matern32.R')
post.dist <- credible_band(t, fet.muf2.m32)
ll <- post.dist$pred - 1.96 * sqrt(diag(post.dist$sigma))
ul <- post.dist$pred + 1.96 * sqrt(diag(post.dist$sigma))

plot(time, muf2(time), type = 'l', lwd = 3, ylim = c(-5,28), xlab = 'Time', ylab = 'Y', cex.lab = 1.5, col = 1)
points(fdata2$x, fdata2$y, pch = 16, cex = .8)
lines(time, post.dist$pred, type = 'l', lwd = 3, col = 5)
lines(time, ll, type = 'l', lwd = 3, col = 5, lty = 2)
lines(time, ul, type = 'l', lwd = 3, col = 5, lty = 2)
legend("topright", c("True", "Matern 3/2"), col = c(1,5), lty = 1, bty = 'n', horiz = F, lwd = 3)
dev.off()