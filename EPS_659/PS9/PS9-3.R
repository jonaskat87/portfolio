# part (a)
N <- 1000
xdata <- rnorm(N,mean=0,sd=1)
spec <- fft(xdata,inverse=TRUE)/length(xdata)
perid <- abs(spec)^2
fpad <- (0:(N-1))/N
plot(fpad,perid,pch=20,col='lightblue',
     main='Periodogram of white noise with zero mean and unit variance',
     xlab='f_m',ylab='|Y_m|^2')

# sample periodogram many times to get empirical average
n <- 5000
data <- matrix(1:(N*n),nrow=N)
for (i in 1:n) {
  xdata <- rnorm(N,mean=0,sd=1)
  spec <- fft(xdata,inverse=TRUE)/length(xdata)
  data[,i] <- abs(spec)^2
}
plot(fpad,rowMeans(data),pch=20,col='lightblue',
     main='Periodogram of white noise with zero mean and unit variance,
     averaged over n samples',xlab='f_m',ylab='<|Y_m|^2>')

# part (b)
xdata <- rnorm(N,mean=0,sd=1)
spec <- fft(xdata,inverse=TRUE)/length(xdata)
perid <- abs(spec)^2 

# plot chi-squared PDF on histogram of time series
r <- hist(perid,xlim=c(0.0,max(perid)),breaks=100,freq=TRUE,
          xlab='<|Y_m|^2>',ylab='Frequency (out of N)',
          main='Histogram of periodogram values')
# get chi-squared PDF with expected mean 1/N and
# rescaled for the histogram.
xx <- seq(0,10,by=0.01)
# rescale histogram by minimizing distance between maxima of cells
# and rescaled PDF values evaluated at midpoints of said cells
d <- dchisq(2*N*(r$mids),df=2)
factor <- sum(d*(r$counts))/sum(d*d)
yy <- factor*dchisq(xx,df=2)
lines(xx/(2*N),yy)

# part (f)
val <- -2*log(1-c(0.9,0.95,0.99,0.999))
val <- sqrt(2*val/N)
A90 <- val[1]
A95 <- val[2]
A99 <- val[3]
A99.9 <- val[4]
f90 <- 200 # frequency associated with amplitude A90
f95 <- 25 # frequency associated with amplitude A95
f99 <- 40 # frequency associated with amplitude A99
f99.9 <- 125 # frequency associated with amplitude A99.9
s90 <- A90*sin(2*pi*f90*fpad) # signal associated with A90
s95 <- A95*cos(2*pi*f95*fpad) # signal associated with A95
s99 <- A99*sin(2*pi*f99*fpad) # signal associated with A99
s99.9 <- A99.9*cos(2*pi*f99.9*fpad) # signal associated with A99.9
data <- xdata+s90+s95+s99+s99.9 # full deterministic signal
specnew <- fft(data,inverse=TRUE)/length(data)
peridnew <- abs(specnew)^2 
# plot new periodogram
plot(fpad,peridnew,pch=20,col='lightblue',
     main='Periodogram of Gaussian noise with deterministic signal',
     xlab='f_m',ylab='|Y_m|^2')
# get percentiles corresponding to aforementioned frequencies
sortperidnew <- sort(peridnew)
percentiles <- rep(0,4)
amplitudes <- (c(A90,A95,A99,A99.9)^2)/4
for (i in 1:4) {
  dist <- (sortperidnew-amplitudes[i])^2
  percentiles[i] <- which.min(dist)
} 
percentiles <- percentiles/N
cat("Percentiles:", percentiles)

# part (g)
perid3 <- (perid[1:(N-2)]+perid[2:(N-1)]+perid[3:N])/3
plot(fpad[2:(N-1)],perid3,pch=20,col='deepskyblue',
     main='Smoothed periodogram of white noise',
     xlab='f_m',ylab='|Y_m|^2')
peridnew3 <- (peridnew[1:(N-2)]+peridnew[2:(N-1)]+peridnew[3:N])/3
plot(fpad[2:(N-1)],peridnew3,pch=20,col='deepskyblue',
     main='Smoothed periodogram of white noise
     with deterministic signal',
     xlab='f_m',ylab='|Y_m|^2')

# part (h)
r <- hist(perid3,xlim=c(0.0,max(perid3)),breaks=100,freq=TRUE,
          xlab='<|Y_m|^2>',ylab='Frequency (out of N)',
          main='Histogram of smoothed periodogram values')
# get chi-squared PDF with expected mean 1/N and
# rescaled for the histogram.
xx <- seq(0,20,by=0.01)
# rescale histogram by minimizing distance between maxima of cells
# and rescaled PDF values evaluated at midpoints of said cells
d <- dchisq(6*N*(r$mids),df=6)
factor <- sum(d*(r$counts))/sum(d*d)
yy <- factor*dchisq(xx,df=6)
lines(xx/(6*N),yy)