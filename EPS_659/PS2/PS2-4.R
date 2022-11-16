# change the following each time:
N = 100
p = 0.1
xq <- 0:N

# plot the binomial PDF
pdf = dbinom(xq, N, p, log = FALSE)
plot(xq, pdf, type="p", main="Binomial PDF",
     sub=paste("N =",N,"and p =",p),
     xlab="k", ylab="p(k)",
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
# plot the binomial CDF
cdf = pbinom(xq, N, p, log = FALSE)
plot(xq, cdf, type="p", main="Binomial CDF",
     sub=paste("N =",N,"and p =",p),
     xlab="k", ylab="p(k)",
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)