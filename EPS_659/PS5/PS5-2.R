# part (a)
Fq <- (0:(10*1000)) / 1000 # sample values of F to plot PDF
# values for degrees of freedom
dfq <- matrix(c(2,8,10,20,17,2000), nrow=3, byrow=TRUE)
for (i in 1:3) {
  Fpdf <- df(Fq, dfq[i,1], dfq[i,2], ncp=0, log=FALSE)
  # linear-linear plot for F PDF
  par(mar = c(5, 5, 5, 5))
  plot(Fq,Fpdf,
       main=paste("F PDF with df1=", dfq[i,1],
                  "and df2=", dfq[i,2]),
       sub="Linear-linear Plot",
       xlab="f", ylab="p(F=f)", col="gold", type="l", 
       cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
  grid()
  # linear-log plot for F PDF
  plot(Fq,Fpdf,log="y",
       main=paste("F PDF with df1=", dfq[i,1],
                  "and df2=", dfq[i,2]),
       sub="Linear-log Plot",
       xlab="f", ylab="p(F=f)", col="mediumblue", type="l", 
       cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
  grid()
}

# part (b)
# load data
fileglobal <- paste0("C:/Users/jonas/Documents",
                    "/EPS 659/PS5/Hadcrut_GlobalAverage.csv")
had_ns_avg <- read.csv(fileglobal)
temps <- had_ns_avg$DTEMPC
years <- had_ns_avg$YEAR
months <- had_ns_avg$MONTH

grandmean <- mean(temps) # grand mean
# index vector for indices at beginning of each decade
indices <- seq(0,length(temps)-12*10,by=12*10)
ms <- rep(0,length(indices)) # decadal means
vs <- rep(0,length(indices)) # decadal variances
Yb <- 0
Yw <- 0
# loop through adjacent decades
for (i in 1:length(indices)) {
  ms[i] <- mean(temps[1:(12*10)+indices[i]]) # decadal mean
  vs[i] <- var(temps[1:(12*10)+indices[i]]) # decadal variance
  Yb <- Yb + (ms[i]-grandmean)^2 # compute sum for Yb
  # compute outer sum for Yw (inner sum already computed)
  Yw <- Yw + vs[i] 
}
Yb <- Yb / (length(indices)-1) # scale Yb by 1/(N-1)
Yw <- Yw / length(indices) # scale Yw by 1/N
Fratio <- Yb / Yw # compute F (variance ratio)

stds <- sqrt(vs) # standard deviations of the sample-means
decmids <- 10*(0:(length(indices)-1)) + 1855 # vector of decadal midpoints
# plot decadal means and variances
par(mar = c(5, 5, 5, 5))
plot(decmids, ms, type="l", col="red", 
     ylim=range(c(ms-stds, ms+stds)),
     main="Sample Mean and Variance for Temperature Anomalies,
     Decade-by-decade",
     xlab="t [years]", ylab="Value of Statistic",
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
lines(decmids, vs, type="l", col="blue")
# add error bars
arrows(x0=decmids, y0=ms-stds, x1=decmids, y1=ms+stds,
       code=3, angle=90, length=0.05, col="dimgrey", lwd=1)
grid()

print(paste("Yb (decadal):", Yb))
print(paste("Yw (decadal):", Yw))
print(paste("F (variance ratio, decadal):", Fratio))
print(paste("Statistical probability for nonrandomness (decadal):", 
            pf(Fratio, length(indices)-1,
               length(indices)*(12*10-1))))
print(paste("Complementary probability of a Type 1 error (decadal):", 
            1-pf(Fratio, length(indices)-1,
                 length(indices)*(12*10-1))))

# part (c)
# monthly temperature data for 1850-1970
# sorted into a matrix with 12 columns
monthtemp <- matrix(temps[1:(12*10*12)], ncol=12, byrow=TRUE)
grandmean <- mean(temps[1:(12*10*12)]) # grand mean for 1850-1970
ms <- rep(0,12) # monthly means
vs <- rep(0,12) # monthly variances
Yb <- 0
Yw <- 0
# loop through months
for (i in 1:12) {
  ms[i] <- mean(monthtemp[,i]) # monthly mean
  vs[i] <- var(monthtemp[,i]) # monthly variance
  Yb <- Yb + (ms[i]-grandmean)^2 # compute sum for Yb
  # compute outer sum for Yw (inner sum already computed)
  Yw <- Yw + vs[i] 
}
Yb <- Yb / (12-1) # scale Yb by 1/(N-1)
Yw <- Yw / 12 # scale Yw by 1/N
Fratio <- Yb / Yw # compute F (variance ratio)

stds <- sqrt(vs) # standard deviations of the sample-means
# plot monthly means and variances
par(mar = c(5, 5, 5, 5))
plot(1:12-0.5, ms, type="l", col="springgreen3", 
     ylim=range(c(ms-stds, ms+stds, vs)),
     main="Sample Mean and Variance for Temperature Anomalies,
     Month-by-month",
     xlab="t [months]", ylab="Value of Statistic",
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
lines(1:12-0.5, vs, type="l", col="lightsalmon2")
# add error bars
arrows(x0=1:12-0.5, y0=ms-stds, x1=1:12-0.5, y1=ms+stds,
       code=3, angle=90, length=0.05, col="dimgrey", lwd=1)
grid()

print(paste("Yb (monthly):", Yb))
print(paste("Yw (monthly):", Yw))
print(paste("F (variance ratio, monthly):", Fratio))
print(paste("Statistical probability for nonrandomness (monthly):", 
            pf(Fratio,  12-1,
               12*(120-1))))
print(paste("Complementary probability of a Type 1 error (monthly):", 
            1-pf(Fratio,  12-1,
                 12*(120-1))))
