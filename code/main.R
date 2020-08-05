rm(list=ls())
library(readr)
library(MASS)
library(kernlab)
library(mvtnorm)
library(Matrix)
library(optimx)

source("code/functions.R")


bone <- read_csv("data/spnbmd.csv")
bone <- as.data.frame(bone)


# data without discarding subjects with single observation

q <- table(bone$idnum)
q <- row.names(q)[as.numeric(table(bone$idnum))>1]
bone <- bone[bone$idnum %in% q,]


# Split data by class labels

train <- split(bone, bone$sex) # for classifying sex

plot(1, type="n", xlab = 'Time', ylab = 'Spinal Bone Mineral Density', 
     xlim = c(min(bone$age), max(bone$age)),
     ylim = c(min(bone$spnbmd), max(bone$spnbmd)),
     cex.main = 1.5, cex.lab = 1.5)
for(i in bone$idnum){
  Cr <- ifelse(bone[bone$idnum == i,4] == 'mal', 1, 0)
  x <- as.numeric(bone[bone$idnum == i,3])
  y <- as.numeric(bone[bone$idnum == i,5])
  lines(x, y, type = 'b', 
        lwd = 4, col = rgb(Cr, 0, 0, .3))
}
legend('bottomright', c('Male', 'Female'), lty = 1, lwd = 2, bty = 'n', 
       col = c(rgb(1, 0, 0, 1),col = rgb(0, 0, 0, 1)))










#choose kernel
source("code/RBF.R")
source("code/laplace.R")
source("code/matern52.R")
source("code/matern32.R")



# Fit the model
gpFit.male <- feature(id = train$mal$idnum, x = train$mal$age, y = train$mal$spnbmd)
gpFit.female <- feature(id = train$fem$idnum, x = train$fem$age, y = train$fem$spnbmd)

save(list = c('gpFit.male', 'gpFit.female'), file = './data/bone/laplace.Rdata')

testTime <- seq(min(bone$age), max(bone$age), length.out = 100)

postDist.male <- postDist(testTime, gpFit.male)
postDist.female <- postDist(testTime, gpFit.female)


# plot
plot(1, type="n", xlab = 'Age (Years)', ylab = 'Spinal Bone Mineral Density', 
     xlim = c(min(bone$age), max(bone$age)),
     ylim = c(min(bone$spnbmd), max(bone$spnbmd)),
     cex.main = 1.5, cex.lab = 1.5)
for(i in bone$idnum){
  Cr <- ifelse(bone[bone$idnum == i,4] == 'mal', 1, 0)
  x <- as.numeric(bone[bone$idnum == i,3])
  y <- as.numeric(bone[bone$idnum == i,5])
  lines(x, y, type = 'b', 
        lwd = 1, col = rgb(Cr, 0, 0, .3))
}
lines(testTime, postDist.male$mu, lwd = 4, col = 2)
lines(testTime, postDist.female$mu, lwd = 4, col = 1)

lines(testTime, postDist.male$ul, lwd = 3, lty = 2, col = 2)
lines(testTime, postDist.male$ll, lwd = 3, lty = 2, col = 2)

lines(testTime, postDist.female$ul, lwd = 3, lty = 2, col = 1)
lines(testTime, postDist.female$ll, lwd = 3, lty = 2, col = 1)

legend('bottomright', c('Male', 'Female'), lty = 1, lwd = 2, bty = 'n', 
       col = c(rgb(1, 0, 0, 1),col = rgb(0, 0, 0, 1)))






