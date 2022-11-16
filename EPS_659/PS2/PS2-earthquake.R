centuries <- 100 # number of centuries to realize
N <- 100 # number of years per century
p <- 0.01 # earthquake probability/year
M <- 100000 # number of samples to take 

data <- rep(0, M)
# sample 100-century realizations
for (i in 1:M) {
  data[i] = sum(rbinom(centuries, N, p) >= 3, na.rm=TRUE)
}
tmean <- N*0.0793732
# plot histogram of results
# how many centuries with K>=3 earthquakes were observed
# during each 100-century period?
hist(data, main="Earthquake cluster frequency",
     xlab="Number of centuries with K>=3 earthquakes",
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
abline(v=tmean, col="purple",lwd=1)

smean = mean(data)
print(smean)
print(tmean)
