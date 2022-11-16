# Start the clock!
ptm <- proc.time()

k <- 100 # set number of cases
n <- 1000 # set number of realizations for each case
mx <- 1 # mean for X
sx <- 2 # standard deviation for X
my <- -1 # mean for Y
sy <- 0.5 # standard deviation for Y

cases <- seq(1, k, by=1)
stat <- matrix(0, nrow=6, ncol=k)
exact <- rep(0, 6)
samplem <- rep(0, 6)
samplev <- rep(0, 6)
exact[1] <- mx
exact[2] <- mx^2 + sx^2
exact[3] <- my
exact[4] <- my^2 + sy^2
exact[5] <- 0.5 * (mx + my)
exact[6] <- 0.25 * (mx^2 + my^2 + sx^2 + sy^2) + 
  0.5*mx*my
for (i in cases) {
  x <- rnorm(n, mx, sx)
  y <- rnorm(n, my, sy)
  z <- (x+y)/2
  stat[1, i] <- mean(x)
  stat[2, i] <- mean(x^2)
  stat[3, i] <- mean(y)
  stat[4, i] <- mean(y^2)
  stat[5, i] <- mean(z)
  stat[6, i] <- mean(z^2)
}
for (j in seq(1, 6, by=1)) {
  hist(stat[j,], main=paste("Statistic",j), 
       xlab="Sample value", breaks=10,
       cex.lab=2.0, cex.axis=2.0, cex.main=2.0, cex.sub=2.0)
  abline(v=exact[j],col="blue",lwd=2)
  samplem[j] <- mean(stat[j,])
  samplev[j] <- var(stat[j,])
}
print(samplem)
print(samplev)

# Stop the clock
print(proc.time() - ptm)