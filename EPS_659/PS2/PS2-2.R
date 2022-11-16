# part b):
nu <- 4 # degrees of freedom
a <- 0
b <- 10
dx <- 0.001 # step-size (for integral too)
xq <- seq(a, b, by=dx) # x-values

pq = dchisq(xq, nu) # PDF
plot(xq, pq, type="l", main="PDF for the chi-squared distribution", 
     sub=paste("degrees of freedom =",nu),
     xlab="x", ylab="p(x)", cex.lab=1.5, 
     cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
abline(v=nu-2, col="green", lwd=2)

# part c):
ub <- 50 # upper bound for numerical integral
sample <- seq(0, ub, by=dx)
# initialize matrix of integrand values
intq <- matrix(0, nrow=2, ncol=length(sample)) 
intq[1,] <- dchisq(sample, nu)
intq[2,] <- sample*intq[1,]
# composite trapezoidal rule simultaneously for both integrals
int <- dx*intq[,1]/2
for (i in 2:(length(sample)-1)) {
  int <- int + dx*intq[,i]
}
int <- int + (dx*intq[,length(sample)]/2)

# part d):
nu <- 10 # degrees of freedom
b <- 30
xq <- seq(a, b, by=dx) 

pq = dchisq(xq, nu)
plot(xq, pq, type="l", main="PDF for chi-squared distribution", 
     sub=paste("degrees of freedom =",nu),
     xlab="x", ylab="p(x)", cex.lab=1.5, 
     cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
abline(v=nu-2, col="green", lwd=2)

# part e):
nuq <- seq(2, 10, 2)
Pq <- rep(0, 5)
for (i in 1:5) {
  Pq[i] <- pchisq(2*nuq[i], nuq[i])
}
Pq <- 1 - Pq # we want P(X>2*nu), not P(X<2*nu)

# results
print(ub/dx) # print number of elements for integral
print(int)
print(Pq)
