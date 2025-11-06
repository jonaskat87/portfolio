dx <- 1e-4 # sub-interval size
tol <- 1e-10 # relative error tolerance
exact <- 1+log(2)+digamma(1) # exact value for integral
maxM <- 1e6 # max number of points (to prevent overflow)

intq <- rep(0, maxM) # value of integral as M increases
errq <- rep(1, maxM) # value of relative error as M increases

# composite trapezoidal rule
intq[2] <- dx*(0.25*dx*log(dx)*exp(-dx/2))/2
errq[2] <- abs(intq[2]-exact)/exact
xf <- dx # upper bound for integral
M <- 2 # number of quadrature points
while ((errq[M]>=tol) && (M<maxM)) {
  intq[M+1] <- intq[M]+
    0.5*dx*(0.25*xf*log(xf)*exp(-xf/2)+
              0.25*(xf+dx)*log(xf+dx)*exp(-(xf+dx)/2))
  errq[M+1] <- abs(intq[M+1]-exact)/exact
  M <- M+1
  xf <- xf+dx
}

plot(dx*(0:(M-1)),intq[1:M],type="l",
     main="Convergence of geometric mean integral",
     sub=paste("dx =",dx,"and error tolerance =",1e-10),
     xlab="Upper bound for integration",ylab="Value of integral",
     cex.lab=1.5,cex.axis=1.0,cex.main=1.5,cex.sub=1.5)
plot(dx*(0:(M-1)),errq[1:M],type="l",log='y',
     main="Errors for geometric mean integral (semi-log)",
     sub=paste("dx =",dx,"and error tolerance =",1e-10),
     xlab="Upper bound for integration",ylab="Relative error",
     cex.lab=1.5,cex.axis=1.0,cex.main=1.5,cex.sub=1.5)

print(intq[M])
print(errq[M])
print(xf)
print(M)