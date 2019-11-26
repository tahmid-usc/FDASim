funcgen <- function(muf, z, theta) {
  x <- seq(0,1,length.out = 100)
  mu <- muf(x, z)
  gt <- mu + mvrnorm(1, rep(0,100), ker(x, l = theta[1], sigf = theta[2])) 
  return(list(x=x, gt = gt))
}

muf <- function(x, z) {
  sin(x*10 + z)
}


func <- funcgen(muf, z = 0, theta = c(1,.1))
plot(func$x, func$gt, type = 'l', lwd = 3, ylim = c(-2,2))
for(i in 1:10) {
  if(i <= 5){
  func <- funcgen(muf, z = 0, theta = c(1,.1))
  lines(func$x, func$gt, lwd = 3)}
  func <- funcgen(muf, z = 1, theta = c(1,.1))
  lines(func$x, func$gt, lwd = 3, col = 2)
}
legend('bottomright', c('z = 0', 'z = 1'), col = 1:2, bty = 'n', lwd = 2)


fdagen <- function(n = 100, maxt = 20, muf, theta = rep(1,3)) {

  train <- data.frame()
  #n number of functions in the sample data
  n.time <- sample(1:maxt, size=n, replace=T) 
  
  for(i in 1:n) {
    id <- rep(i, n.time[i])
    x <- sort(runif(n.time[i], 0 , 100))
    mu <- muf(x)
    gt <- mvrnorm(1, rep(0, n.time[i]), ker(x, l = theta[1], sigf = theta[2])) 
    y <- mu + gt
    #y <- mu + gt + rnorm(n.time[i], 0, theta[3])
    train <- rbind(train, cbind(x, y, id))
  }
  return(train)
}

muf <- function(x) {
  sin(x/5)
}


fdata <- fdagen(n = 50, maxt = 10, muf = muf, theta = c(.1,.2,.01))

#plotting FDA
t <- 0:100
plot(t, muf(t), type = 'l', lwd = 4, ylim = c(-3,3))
for(i in 1:100){
  lines(fdata[fdata$id == i,1], fdata[fdata$id == i,2], type = 'b', lwd = 2, col = 'grey')
}



