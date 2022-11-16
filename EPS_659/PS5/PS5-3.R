fileSH <- paste0("C:/Users/jonas/Documents",
                 "/EPS 659/PS5/Hadcrut_SHAverage.csv")
fileNH <- paste0("C:/Users/jonas/Documents",
                 "/EPS 659/PS5/Hadcrut_NHAverage.csv")
had_sh_avg <- read.csv(fileSH)
had_nh_avg <- read.csv(fileNH)
temps_sh <- had_sh_avg$DTEMPC
temps_nh <- had_nh_avg$DTEMPC

# part (a)
# index vector for indices at beginning of each decade
indices <- seq(0,length(temps_nh)-12*10,by=12*10)
Nms <- rep(0,length(indices)) # Northern Hemisphere decadal means
Nvs <- rep(0,length(indices)) # " " " variances
Sms <- rep(0,length(indices)) # Southern Hemisphere decadal means
Svs <- rep(0,length(indices)) # " " " variances
rs <- rep(0,length(indices)) # correlation coefficients
# loop through adjacent decades
for (i in 1:length(indices)) {
  Ntempsi = temps_nh[1:(12*10)+indices[i]]
  Stempsi = temps_sh[1:(12*10)+indices[i]]
  Nms[i] <- mean(Ntempsi) # NH mean
  Nvs[i] <- var(Ntempsi) # NH variance
  Sms[i] <- mean(Stempsi) # SH mean
  Svs[i] <- var(Stempsi) # SH variance
  topfactor = sum((Ntempsi-Nms[i])*(Stempsi-Sms[i]))
  bottomfactor = sum((Ntempsi-Nms[i])^2)*sum((Stempsi-Sms[i])^2)
  rs[i] <- topfactor / sqrt(bottomfactor)
}
Nstds <- sqrt(Nvs)
Sstds <- sqrt(Svs)
# vector of decadal midpoints
decmids <- 10*(0:(length(indices)-1)) + 1855
# plot decadal means with error bars
par(mar = c(5, 5, 5, 5))
plot(decmids, Nms, type="l", col="red", 
     ylim=range(c(Nms-Nstds, Nms+Nstds, Sms-Sstds, Sms+Sstds)),
     main="Decadal Sample Mean for Temperature Anomalies in the
     Northern vs. Southern Hemisphere",
     xlab="Midpoint of Time Interval [years]", 
     ylab="Mean",
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
lines(decmids, Sms, type="l", col="blue")
# add error bars
arrows(x0=decmids, y0=Nms-Nstds, x1=decmids, y1=Nms+Nstds,
       code=3, angle=90, length=0.05, col="red", lwd=1)
arrows(x0=decmids, y0=Sms-Sstds, x1=decmids, y1=Sms+Sstds,
       code=3, angle=90, length=0.05, col="blue", lwd=1)
grid()

# plot correlation coefficients
plot(decmids, rs, type="l", col="palegreen3", 
     main="Decadal Sample Correlations for Temperature Anomalies 
     in the Northern vs. Southern Hemisphere",
     xlab="Midpoint of Time Interval [years]", 
     ylab="Sample Correlation Coefficient",
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
grid()

# part (b)
# confidence for non-randomness for correlations
conf <- pf((12*10-2)*(rs^2), 1, 12*10-2)
# plot confidences
plot(decmids, conf, type="l", col="orchid3", 
     main="Confidence for Non-randomness for
the Correlations in Part (a)",
     xlab="Midpoint of Time Interval [years]", 
     ylab="Confidence Level, P(F<=(N-2)r^2)",
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
grid()

# part (c)
# plot complementary probabilities of a Type 1 error
par(mar = c(7, 7, 7, 7))
plot(decmids, 1-conf, type="l", col="skyblue3", log="y", 
     main="Complementary Probabilities of a Type 1 Error
 for the Correlations in Part (a)",
     xlab="Midpoint of Time Interval [years]", 
     ylab="Probability of a Type 1 error, 
     1-P(F<=(N-2)r^2)=P(F>=(N-2)r^2)",
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
grid()

# part (d)
Ngrandmean <- mean(temps_nh) # NH grand mean
Sgrandmean <- mean(temps_sh) # SH grand mean 
topfactor = sum((Nms-Ngrandmean)*(Sms-Sgrandmean))
bottomfactor = sum((Nms-Ngrandmean)^2)*sum((Sms-Sgrandmean)^2)
rb <- topfactor / sqrt(bottomfactor)
conf <- pf((length(Nms)-1)*rb, 1, length(Nms)-1)
cat(paste("Correlation between the decadal mean temperatures
for each hemispheres:", rb))
cat(paste("Confidence for nonrandomness of correlation 
          between the sample means:", conf))
cat(paste("Complementary probability of a Type 1 error:",
          1-conf))