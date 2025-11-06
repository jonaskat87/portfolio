# load data
file = paste0("C:\\Users\\jonas\\Documents",
"\\EPS 659\\PS3\\sunspots_annual_1700_2008.csv")
sunspots <- read.csv(file)
yr <- sunspots$YEAR # vector of years
sspots <- sunspots$SUNACTIVITY # vector of sunspots

# part a):
hist_info <- hist(sspots, main=
                    paste("Distribution of annual-average sunspots"), 
     xlab="Number of sunspots", ylab="Frequency", breaks=20,
     cex.lab=2.0, cex.axis=1.0, cex.main=2.0, cex.sub=1.0)

# part b):
m = mean(sspots)
s = sqrt(var(sspots))

# part c):
range <- floor(min(sspots)):ceiling(max(sspots))
dd <- dnorm(range,mean=m,sd=s)
# plot PDF and histogram
par(mfrow=c(2,1)) 
plot(range,dd,type='l',main=paste("Gaussian PDF"),
     xlab="s=number of sunspots", ylab="p(s)", 
     cex.lab=2.0, cex.axis=1.0, cex.main=2.0, cex.sub=1.0)
hist(sspots, main=
       paste("Distribution of annual-average sunspots"), 
     xlab="Number of sunspots", ylab="Frequency", breaks=20,
     cex.lab=2.0, cex.axis=1.0, cex.main=2.0, cex.sub=1.0)
# total area underneath Gaussian PDF
PDFarea <- pnorm(ceiling(max(sspots)),mean=m,sd=s)-
  pnorm(floor(min(sspots)),mean=m,sd=s)
# extract area underneath histogram
hist_counts <- hist_info$counts # Extract histogram counts
hist_mids <- hist_info$mids # Extract histogram midpoints of cells
# area under histogram
hist_area <- sum(hist_counts*(hist_mids[2]-hist_mids[1]))
# normalize PDF by area
ddscaled <- dd * hist_area / PDFarea
# plot histogram with scaled PDF
par(mfrow=c(1,1))
hist(sspots, main=paste("Distribution of annual-average sunspots"), 
     xlab="Number of sunspots", ylab="Frequency", breaks=20,
     cex.lab=2.0, cex.axis=1.0, cex.main=2.0, cex.sub=1.0)
lines(range,ddscaled)

# part d):
adjlogs <- log(10+sspots)
hist_infolog <- hist(adjlogs, 
                  main=paste("Distribution of annual-average sunspots"), 
                  xlab="Adjusted logarithmic number of sunspots", 
                  ylab="Frequency", breaks=50,
                  cex.lab=2.0, cex.axis=1.0, cex.main=2.0, cex.sub=1.0)
mlog <- mean(adjlogs)
slog <- sqrt(var(adjlogs))

# part e):
rangelog <- seq(floor(min(adjlogs)),ceiling(max(adjlogs)),by=0.01)
ddlog <- dnorm(rangelog,mean=mlog,sd=slog)
# plot PDF and histogram
par(mfrow=c(2,1)) 
plot(rangelog,ddlog,type='l',main=paste("Gaussian PDF (log S)"),
     xlab="log(10+s)", ylab="p(log(10+s))", 
     cex.lab=2.0, cex.axis=1.0, cex.main=2.0, cex.sub=1.0)
hist(adjlogs, main=
       paste("Distribution of annual-average sunspots"), 
     xlab="Adjusted logarithmic number of sunspots", 
     ylab="Frequency", breaks=50,
     cex.lab=2.0, cex.axis=1.0, cex.main=2.0, cex.sub=1.0)
# total area underneath Gaussian PDF
PDFarealog <- pnorm(ceiling(max(adjlogs)),mean=mlog,sd=slog)-
  pnorm(floor(min(adjlogs)),mean=mlog,sd=slog)
# extract area underneath histogram
hist_countslog <- hist_infolog$counts # Extract histogram counts
hist_midslog <- hist_infolog$mids # Extract cell midpoints
# area under histogram
hist_arealog <- sum(hist_countslog*(hist_midslog[2]-hist_midslog[1]))
# normalize PDF by area
ddscaledlog <- ddlog * hist_arealog / PDFarealog
# plot histogram with scaled PDF
par(mfrow=c(1,1))
hist(adjlogs, 
     main=paste("Distribution of annual-average sunspots"), 
     xlab="Adjusted logarithmic number of sunspots", 
     ylab="Frequency", breaks=50,
     cex.lab=2.0, cex.axis=1.0, cex.main=2.0, cex.sub=1.0)
lines(rangelog,ddscaledlog)

# part f):
ssorted <- sort(sspots)
ecdfeval <- ecdf(ssorted)
xq <- seq(floor(min(sspots))-10,ceiling(max(sspots))+10,by=0.01)
plot(xq,ecdfeval(xq),type='l',main=paste("Empirical CDF"),
     xlab="s", ylab="P(S<s)", 
     cex.lab=2.0, cex.axis=1.0, cex.main=2.0, cex.sub=1.0)

# part g):
mus <- rep(0,3)
sigmas <- rep(0,3)
colors <- c("orange", "green", "purple") 
# get mu and sigma from empirical median and empirical mean
mus[1] <- log(median(10+sspots)) 
sigmas[1] <- sqrt(2*(log(mean(10+sspots))-mu1))
# get mu and sigma from empirical mean and standard deviation
mus[2] <- log((mean(10+sspots)^2)/
                sqrt(mean(10+sspots)^2+sd(10+sspots)^2))
sigmas[2] <- sqrt(2)*sqrt(log(sqrt(
  mean(10+sspots)^2+sd(10+sspots)^2)/mean(10+sspots)))
# get mu and sigma by adding a bit to sigma2
mus[3] <- mus[1]
sigmas[3] <- sigmas[2]+0.15
plot(xq+10,ecdfeval(xq),type='l',main=paste("CDF Comparison"),
     xlab="s", ylab="P(S+10<s)",col="red",
     cex.lab=2.0, cex.axis=1.0, cex.main=2.0, cex.sub=1.0)
for (i in 1:3) {
  lines(xq+10,plnorm(xq+10,meanlog=mus[i],sdlog=sigmas[i]),
        type='l',col=colors[i])
}

# outputs for a)
print(max(sspots))
print(yr[which.max(sspots)])
print(min(sspots))
print(yr[which.min(sspots)])

# outputs for b)
print(m)
print(s)

# outputs for d)
print(max(adjlogs))
print(yr[which.max(adjlogs)])
print(min(adjlogs))
print(yr[which.min(adjlogs)])
print(mlog)
print(slog)

# outputs for g)
print(mus)
print(sigmas)