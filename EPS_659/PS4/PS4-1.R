file = paste0("C:/Users/jonas/Documents",
              "/EPS 659/PS4/Hadcrut_GlobalAverage.csv")
had_ns_avg <- read.csv(file)
temps <- had_ns_avg$DTEMPC
years <- had_ns_avg$YEAR
months <- had_ns_avg$MONTH

time <- years + (months-0.5)/12.0
N <- length(time) # number of data points

# part (a)
# computing moving avg
# c(tail(x, -n), rep(NA, n)) shifts x forward by n indices
# c(rep(NA, n), head(x, -n)) shifts x backward by n indices
avg <- temps
# compute sum for running average
for (i in 1:5) {
  avg <- avg + c(tail(temps, -i), rep(NA, i)) +
    c(rep(NA, i), head(temps, -i)) 
}
avg <- avg[6:(N-5)]/11 # filter out undefined values
par(mar = c(5, 5, 5, 5)) # Set the margin on all sides to 5
plot(time, temps, type="l", main="Global-Average Monthly Temperature 
     Anomalies from 1850AD to 2021AD",
     xlab="t [years]", ylab="T [°C]", col="black",
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
lines(time[6:(N-5)], avg, col="lightgreen")
grid()

# part (b)
sm <- mean(temps) # sample mean
sv <- var(temps) # sample variance
print(paste("The sample mean is", sm))
print(paste("The sample variance is", sv))
x <- seq(min(temps), max(temps), by=0.05)
pdf <- dnorm(x, mean=sm, sd=sqrt(sv)) # Gaussian PDF

# plot histogram of data
hist_info <- hist(temps, main="Global-Average Monthly Temperature 
     Anomalies from 1850AD to 2021AD", 
     xlab="T [°C]", ylab="Frequency", breaks=50,
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
# total area underneath Gaussian PDF
PDFarea <- pnorm(max(temps), mean=sm, sd=sqrt(sv))-
  pnorm(min(temps), mean=sm, sd=sqrt(sv))
# extract area underneath histogram
hist_counts <- hist_info$counts # Extract histogram counts
hist_mids <- hist_info$mids # Extract histogram midpoints of cells
# area under histogram
hist_area <- sum(hist_counts*(hist_mids[2]-hist_mids[1]))
# normalize PDF by area
pdfscaled <- pdf * hist_area / PDFarea
# also plot PDF on histogram, scaled as suggested in problem
lines(x, pdfscaled, type="l", col="orangered")
grid()

# part (c)
tempsorted <- sort(temps) # sort temp data
ecdfeval <- ecdf(tempsorted) # empirical CDF of temp data
# evaluate empirical CDF at sample points
empiricalCDF <- ecdfeval(tempsorted) 
# " Gaussian CDF """
gaussianCDF <- pnorm(tempsorted, mean=sm, sd=sqrt(sv)) 
# Q-Q plot
par(mar = c(7, 7, 7, 7))
plot(gaussianCDF, empiricalCDF, type="l", col="turquoise", 
    main="Q-Q Plot for Empirical CDF vs. Gaussian CDF",
     xlab="Normal Theoretical Quantiles", ylab=
      paste0("Empirical Quantiles", "for Hadley-Center 
     Global Average Temperature Anomalies"),
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
lines(x, x, type="l", lty=2, col="black")
grid()
# find largest deviation
nonsortedGaussianCDF <- pnorm(temps, mean=sm, sd=sqrt(sv))
devs <- abs(ecdfeval(temps)-
             nonsortedGaussianCDF) 
print(max(devs)) # what is the largest deviation?
print(temps[which.max(devs)]) # at which temp value?
print(nonsortedGaussianCDF[which.max(devs)]) # at which quantile?
print(time[which.max(devs)]) # at which year?

# part (d)
med <- median(tempsorted) # median temperature
# compute rho (difference between midpoints of empirical and Gaussian CDFs)
rho <- pnorm(med, mean=sm, sd=sqrt(sv)) - 0.5 
betafactor <- (1-4*rho^2)^(N/2) # compute factor
print(betafactor) # return factor

# part (e)
# there are 12*10 months in a decade
# there are floor(N/(12*10)) full decades in the data
factor <- 12*10
count <- floor(N/factor)
dms <- rep(0, count) # sample means for each decade
dvs <- rep(0, count) # sample variances """
for (i in 1:count) {
  dms[i] <- mean(temps[(factor*(i-1)+1):(factor*i)])
  dvs[i] <- var(temps[(factor*(i-1)+1):(factor*i)])
}
dstds <- sqrt(dvs) # standard deviations of the sample-means
decmids <- 10*(0:(count-1)) + 1855 # vector of decadal midpoints
# plot decadal means and variances
par(mar = c(5, 5, 5, 5))
plot(decmids, dms, type="l", col="red", 
     ylim=range(c(dms-dstds, dms+dstds)),
     main="Sample Mean and Variance for Temperature Anomalies,
     Decade-by-decade",
     xlab="t [years]", ylab="Value of Statistic",
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
lines(decmids, dvs, type="l", col="blue")
# add error bars
arrows(x0=decmids, y0=dms-dstds, x1=decmids, y1=dms+dstds,
       code=3, angle=90, length=0.05, col="dimgrey", lwd=1)
grid()

# part (f)
# by solving 10(x-1)+1855=y for 1905 and 1995, we find that the 
# first decade of the 20th century corresponds to x=6 and the
# last decade of the 20th century corresponds to x=15.
twoms <- c(dms[6], dms[15]) # means for first and last
twovs <- c(dvs[6], dvs[15]) # variances """"
# get indices for full time series data in each decade
indices <- matrix(c((factor*(6-1)+1):(factor*6),
                    (factor*(15-1)+1):(factor*15)), 
                  nrow=2, ncol=factor, byrow = TRUE) 
# sort temperatures in each decade
twotempsort <- matrix(c(sort(temps[indices[1,]]), 
                        sort(temps[indices[2,]])), 
                      nrow=2, ncol=factor, byrow = TRUE) 
labels <- c('First', 'Last') # labels for each statistic
# preallocate Gaussian PDF data
twogPDFs <- matrix(0, nrow=2, ncol=factor) 
# preallocate Gaussian CDF data
twogCDFs <- matrix(0, nrow=2, ncol=factor)
# preallocate empirical CDF data
twoeCDFs <- matrix(0, nrow=2, ncol=factor)
for (i in 1:2) {
  # generate Gaussian PDF 
  twogPDFs[i,] <- dnorm(twotempsort[i,],
                        mean=twoms[i], sd=sqrt(twovs[i]))
  # generate Gaussian CDF
  twogCDFs[i,] <- pnorm(twotempsort[i,],
                        mean=twoms[i], sd=sqrt(twovs[i]))
  # generate empirical CDF 
  ecdfeval <- ecdf(twotempsort[i,])
  twoeCDFs[i,] <- ecdfeval(twotempsort[i,]) 
}

# plot PDFs against each other
plot(twotempsort[1,], twogPDFs[1,], type="l", col="blue", 
     xlim=range(c(twotempsort[1,], twotempsort[2,])),
     main="Gaussian PDF for First Decade (blue) vs. Last Decade (gold)",
     xlab="T [°C]", ylab="p(t)",
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
lines(twotempsort[2,], twogPDFs[2,], type="l", col="gold")
abline(v=twoms[1],col="dimgrey",lwd=1, lty=2)
abline(v=twoms[2],col="dimgrey",lwd=1, lty=2)
grid()

# plot CDFs for first and second decades
colors <- c("blue", "gold")
for (i in 1:2) {
  plot(twogCDFs[i,], twoeCDFs[i,], type="l", col=paste0(colors[i]), 
       main=paste("Q-Q Plot for Empirical CDF vs. Gaussian CDF
       in the ", labels[i], " Decade"),
       xlab="Normal Theoretical Quantiles", ylab="Empirical Quantiles",
       cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
  lines(x, x, type="l", lty=2, col="black")
  grid()
  med <- median(twotempsort[i,]) # median temperature
  # compute rho (difference between midpoints of empirical and Gaussian CDFs)
  rho <- pnorm(med, mean=twoms[i], sd=sqrt(twovs[i])) - 0.5 
  betafactor <- (1-4*rho^2)^(N/2) # compute factor
  print(paste("For the", labels[i], "case:", betafactor)) # return factor
}

# plot theoretical and empirical CDFs against each other
plot(twotempsort[1,], twogCDFs[1,], type="l", col="blue", 
     xlim=range(c(twotempsort[1,], twotempsort[2,])),
     main="Gaussian (solid) and Empirical (dashed) CDFs
     for First Decade (blue) vs. Last Decade (gold)",
     xlab="T [°C]", ylab="P(T<t)",
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
lines(twotempsort[1,], twoeCDFs[1,], type="l", col="blue", lty=2)
lines(twotempsort[2,], twogCDFs[2,], type="l", col="gold")
lines(twotempsort[2,], twoeCDFs[2,], type="l", col="gold", lty=2)
grid()
L1error1 = sum(abs(twoeCDFs[1,]-twogCDFs[1,])) / sum(abs(twogCDFs[1,]))
print(paste("The relative L1 error for the first decade is", L1error1))
L2error1 = sum((twoeCDFs[1,]-twogCDFs[1,])^2) / sum((twogCDFs[1,])^2)
print(paste("The relative L2 error for the first decade is", L2error1))
L1error2 = sum(abs(twoeCDFs[2,]-twogCDFs[2,])) / sum(abs(twogCDFs[2,]))
print(paste("The relative L1 error for the first decade is", L1error2))
L2error2 = sum((twoeCDFs[2,]-twogCDFs[2,])^2) / sum((twogCDFs[2,])^2)
print(paste("The relative L2 error for the first decade is", L2error2))
  
print(paste("Mean for the first decade:", twoms[1]))
print(paste("Mean for the last decade:", twoms[2]))
print(paste("Variance for the first decade:", twovs[1]))
print(paste("Variance for the last decade:", twovs[2]))

