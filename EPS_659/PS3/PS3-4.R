# part a):
k <- 4
s <- 0.5
Nq <- 4:20
nbinq <- (Nq-1)*(Nq-2)*(Nq-3)/(6*2^Nq)
print(nbinq)
plot(Nq,nbinq,
     main="Negative Binomial Distribution",
     sub=paste("k (number of times) =",k,"and s (probability) =",s),
     xlab="N (number of trials)",
     ylab="Probability of kth occurence in the Nth trial",
     cex.lab=1.5,cex.axis=1.0,cex.main=1.5,cex.sub=1.5)

# part b):
k <- 4
s <- 0.5
Nq <- 4:20
nbinqR <- dnbinom(Nq-k, size=k, prob=s)
print(nbinqR)
cumnbinq <- pnbinom(Nq-k, size=k, prob=s)
print(cumnbinq)
Ntoinf <- k:100
plot(Ntoinf,pnbinom(Ntoinf-k, size=k, prob=s),
     main="Negative Binomial Distribution",
     sub=paste("k (number of times) =",k,"and s (probability) =",s),
     xlab="N (number of trials)",
     ylab="Probability of kth occurence in the Nth trial",
     cex.lab=1.5,cex.axis=1.0,cex.main=1.5,cex.sub=1.5)
last10 <- pnbinom(90:100, size=k, prob=s)
print(last10)

# part c):
k <- 4
s <- 0.1
Nq <- k:100
nbinqR <- dnbinom(Nq-k, size=k, prob=s)
print(nbinqR)
cumnbinq <- pnbinom(Nq-k, size=k, prob=s)
print(cumnbinq)
Ntoinf <- k:300
plot(Ntoinf,pnbinom(Ntoinf-k, size=k, prob=s),
     main="Negative Binomial Distribution",
     sub=paste("k (number of times) =",k,"and s (probability) =",s),
     xlab="N (number of trials)",
     ylab="Probability of kth occurence in the Nth trial",
     cex.lab=1.5,cex.axis=1.0,cex.main=1.5,cex.sub=1.5)
last10 <- pnbinom(290:300, size=k, prob=s)
print(last10)

