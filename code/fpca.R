dt1 <- MakeFPCAInputs(IDs = train$train1$id, train$train1$x, train$train1$y, sort = FALSE)
dt2 <- MakeFPCAInputs(IDs = train$train2$id, train$train2$x, train$train2$y, sort = FALSE)


pca1 <- FPCA(dt1$Ly, dt1$Lt, optns = list(dataType = 'Sparse'))
pca2 <- FPCA(dt2$Ly, dt2$Lt, optns = list(dataType = 'Sparse'))

plot(pca1)
plot(pca2)



a1 <- 2
b1 <- .1
a2 <- 2.5
b2 <- 1
c2 <- 0.06

mu1 <- a1 + b1 * pca1$workGrid/5
mu2 <- a2 + b2 * sin(pca2$workGrid/5) + c2 * pca2$workGrid 

mean((mu1 - pca1$mu)^2)
mean((mu2 - pca2$mu)^2)

mu1.gp <- gpsmooth(pca1$workGrid, fet$train1)
mu2.gp <- gpsmooth(pca2$workGrid, fet$train2)

#FPCA
mean((mu1 - pca1$mu)^2)
mean((mu2 - pca2$mu)^2)
# GP
mean((mu1 - mu1.gp)^2)
mean((mu2 - mu2.gp)^2)
