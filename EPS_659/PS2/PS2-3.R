# part a):
m <- 0 # mean of log
s <- 1 # standard deviation of log
dx <- 0.01 # sub-interval size
xq <- seq(0, 6, by=dx)

# plot the log-normal PDF
pdf = dlnorm(xq, m, s, log = FALSE)
plot(xq, pdf, type="l", main="Log-normal PDF",
     sub=paste("mu =",m,"and sigma =",s),
     xlab="x", ylab="p(x)",
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
abline(v=exp(m + 0.5*s^2),col="blue",lwd=1)
# plot the log-normal CDF
cdf = plnorm(xq, m, s, log = FALSE)
plot(xq, cdf, type="l", main="Log-normal CDF",
     sub=paste("mu =",m,"and sigma =",s),
     xlab="x", ylab="P(X<x)",
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
abline(v=exp(m + 0.5*s^2),col="blue",lwd=1)

# part b):
ub <- 500 # upper bound for numerical integral
sample <- seq(0, ub, by=dx)
# initialize matrix of integrand values
intq <- matrix(0, nrow=2, ncol=length(sample)) 
intq[1,] <- dlnorm(sample, m, s, log = FALSE)
intq[2,] <- sample*intq[1,]
# composite trapezoidal rule simultaneously for both integrals
intb <- dx*intq[,1]/2
for (i in 2:(length(sample)-1)) {
  intb <- intb + dx*intq[,i]
}
intb <- intb + (dx*intq[,length(sample)]/2)

# part c):
m <- 2 # mean of log
s <- 1 # standard deviation of log
dx <- 0.01 # sub-interval size
xq <- seq(0, 25, by=dx)

# plot the log-normal PDF
pdf = dlnorm(xq, m, s, log = FALSE)
plot(xq, pdf, type="l", main="Log-normal PDF",
     sub=paste("mu =",m,"and sigma =",s),
     xlab="x", ylab="p(x)",
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
abline(v=exp(m + 0.5*s^2),col="blue",lwd=1)
# plot the log-normal CDF
cdf = plnorm(xq, m, s, log = FALSE)
plot(xq, cdf, type="l", main="Log-normal CDF",
     sub=paste("mu =",m,"and sigma =",s),
     xlab="x", ylab="P(X<x)",
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
abline(v=exp(m + 0.5*s^2),col="blue",lwd=1)

# part d):
ud <- 2000 # upper bound for numerical integral
sample <- seq(0, ud, by=dx)
# initialize matrix of integrand values
intq <- matrix(0, nrow=2, ncol=length(sample)) 
intq[1,] <- dlnorm(sample, m, s, log = FALSE)
intq[2,] <- sample*intq[1,]
# composite trapezoidal rule simultaneously for both integrals
intd <- dx*intq[,1]/2
for (i in 2:(length(sample)-1)) {
  intd <- intd + dx*intq[,i]
}
intd <- intd + (dx*intq[,length(sample)]/2)

# results
# from part (b)
print(ub/dx) # print number of elements for integral
print(intb)
# from part (d)
print(ud/dx) # print number of elements for integral
print(intd)