m <- 1 # mean
s <- 2 # standard deviation
dx <- 0.01 # sub-interval size

# part a):
xq <- seq(-6, 6, by=dx)

# plot the Gaussian PDF
pdf = dnorm(xq, m, s)
plot(xq, pdf, type="l", main="Gaussian PDF",
     sub=paste("mu =",m,"and sigma =",s),
     xlab="x", ylab="p(x)",
     cex.lab=2.0, cex.axis=2.0, cex.main=2.0, cex.sub=2.0)
abline(v=m,col="red",lwd=1)
# plot the Gaussian CDF
cdf = pnorm(xq, m, s)
plot(xq, cdf, type="l", main="Gaussian CDF",
     sub=paste("mu =",m,"and sigma =",s),
     xlab="x", ylab="P(X<x)",
     cex.lab=2.0, cex.axis=2.0, cex.main=2.0, cex.sub=2.0)
abline(v=m,col="red",lwd=1)

# parts b), c), and d):
b <- 100+m # upper bound for numerical integral
a <- -100+m # lower ""
sample <- seq(a, b, by=dx)
# initialize matrix of integrand values
intq <- matrix(0, nrow=3, ncol=length(sample)) 
intq[1,] <- dnorm(sample, m, s)
intq[2,] <- sample*intq[1,]
intq[3,] <- ((sample - m)^2)*intq[1,]
# composite trapezoidal rule simultaneously for all three integrals
int <- dx*intq[,1]/2
for (i in 2:(length(sample)-1)) {
  int <- int + dx*intq[,i]
}
int <- int + (dx*intq[,length(sample)]/2)

# parts f):
lb <- 10
ub <- 90

print((b - a)/dx) # print number of elements
print(int) # print numerical integral
print(qnorm(0.01*c(lb, ub), m, s))