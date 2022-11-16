# part (a)
file <- paste0("C:/Users/jonas/Documents",
               "/EPS 659/PS9/co2_maunaloa.csv")
MaunaLoa <- read.csv(file)
year <- MaunaLoa$YEAR
month <- MaunaLoa$MONTH
time <- MaunaLoa$TIME
co2 <- MaunaLoa$CO2

time2 <- time^2
lm2_co2 <- lm(co2~time+time2)
s2 <- summary(lm2_co2)
y2 <- s2$coefficients[3,1] * time2 + 
  s2$coefficients[2,1] * time + s2$coefficients[1,1]
co2_resid <- co2-y2

# plot quadratic model and raw data vs. time
par(mar=c(5,5,5,5))
plot(time,co2,"l",main="Mauna Loa Carbon Dioxide (ppm),
     Quadratic Model",
     col="cyan3", xlab="Time [years]", 
     ylab="CO2 Per Unit Air [ppm]",
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
lines(time,y2, col="blue")

# plot residuals for quadratic model
plot(time,co2_resid,main="Residuals for Quadratic Model (ppm)",
     col="magenta1", xlab="Time [years]", 
     ylab="Data Minus Model [ppm]",
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5) 

# part (b)
N <- length(co2_resid)
spec <- fft(co2_resid,inverse=TRUE)/N
pgram <- abs(spec)^2
ymax <- max(pgram)+sort(pgram,decreasing=TRUE)[1]-
  sort(pgram,decreasing=TRUE)[2]
ymin <- sort(pgram)[2]
freq <- 12*(0:(N-1))/N
# non-log plot
plot(freq,pgram,ylim=c(ymin,ymax),xlim=c(-0.05,6.0),
     main='Periodogram of noise in Mauna Loa data since 1958
     (regular scale)',
     xlab='f_m',ylab='|Y_m|^2',col='firebrick')
lines(freq,pgram,col='firebrick')
# log plot
plot(freq,pgram,ylim=c(ymin,ymax),xlim=c(-0.05,6.0),
     main='Periodogram of noise in Mauna Loa data since 1958
     (semi-log scale)',
     xlab='f_m',ylab='|Y_m|^2',log="y",col='red')
lines(freq,pgram,col='red')

# part (c)
pgram3 <- (pgram[1:(N-2)]+pgram[2:(N-1)]+pgram[3:N])/3
plot(freq[2:(N-1)],pgram3,ylim=c(ymin,ymax),xlim=c(-0.05,6.0),
     main='Smoothed periodogram of noise in Mauna Loa data since 1958
     (semi-log scale)',
     xlab='f_m',ylab='|Y_m|^2',log="y",col='orangered')
lines(freq,pgram,col='orangered')

# sample inverse frequency functions
plot(freq[2:(N-1)],pgram3,ylim=c(ymin,ymax),xlim=c(-0.05,6.0),
     main='Smoothed periodogram of noise in Mauna Loa data since 1958
     (semi-log scale)',
     xlab='f_m',ylab='|Y_m|^2',log="y",col='orangered')
lines(freq,pgram,col='orangered')
ff <- freq[10:N]
alpha <- c(0.5,1,1.5,2)
freqfunc <- matrix(rep(0,4*((N/2)-9)),nrow=4)
colors <- c("seagreen","royalblue","purple3","navy")
# for loop to compute best fit for A
for (i in 1:4) {
  freqi <- ff[10:(N/2)-9]^(-alpha[i])
  factor <- sum(freqi*pgram3[10:(N/2)])/sum(freqi*freqi)
  freqfunc[i,] <- factor*freqi
  lines(ff[10:(N/2)-9],freqfunc[i,],col=colors[i])
}

# part (d)
co2_diff <- co2_resid[2:N]-co2_resid[1:(N-1)]
spec <- fft(co2_diff,inverse=TRUE)/(N-1)
pgramdiff <- abs(spec)^2
ymax <- max(pgramdiff)+sort(pgramdiff,decreasing=TRUE)[1]-
  sort(pgramdiff,decreasing=TRUE)[2]
ymin <- sort(pgramdiff)[2]
freqdiff <- 12*(0:(N-2))/(N-1)
# non-log plot
plot(freqdiff,pgramdiff,ylim=c(ymin,ymax),xlim=c(-0.05,6.0),
     main='Periodogram of Mauna Loa noise with first-difference filter
     (regular scale)',
     xlab='f_m',ylab='|Y_m|^2',col='purple')
lines(freqdiff,pgramdiff,col='purple')
# log plot
plot(freqdiff,pgramdiff,ylim=c(ymin,ymax),xlim=c(-0.05,6.0),
     main='Periodogram of Mauna Loa noise with first-difference filter
     (semi-log scale)',
     xlab='f_m',ylab='|Y_m|^2',log="y",col='mediumorchid')
lines(freqdiff,pgramdiff,col='mediumorchid')

pgramdiff3 <- (pgramdiff[1:(N-3)]+pgramdiff[2:(N-2)]+
                 pgramdiff[3:(N-1)])/3
plot(freqdiff[2:(N-2)],pgramdiff3,ylim=c(ymin,ymax),xlim=c(-0.05,6.0),
     main='Smoothed periodogram of Mauna Loa noise 
     before and after first-difference filter',type='l',
     xlab='f_m',ylab='|Y_m|^2',log="y",col='green')
lines(freq[2:(N-1)],pgram3,col='red')

# part (e)
co2pad <- rep(0,6000)
co2pad[1:(N-1)] <- co2_diff
dfpad <- 12.0/6000.0
fpad <- dfpad*(0:5999)
spec <- fft(co2pad,inverse=TRUE)/6000
padpgram <- abs(spec)^2
plot(fpad[(6000*0.7/12):(6000*1.3/12)+1],
     padpgram[(6000*0.7/12):(6000*1.3/12)+1],
     main='Filtered periodogram (first frequency interval)',
     xlab='f_m',ylab='|Y_m|^2',col='darkturquoise')
plot(fpad[(6000*1.7/12):(6000*2.3/12)+1],
     padpgram[(6000*1.7/12):(6000*2.3/12)+1],
     main='Filtered periodogram (second frequency interval)',
     xlab='f_m',ylab='|Y_m|^2',col='darkorchid')
plot(fpad[(6000*2.7/12):(6000*3.3/12)+1],
     padpgram[(6000*2.7/12):(6000*3.3/12)+1],
     main='Filtered periodogram (third frequency interval)',
     xlab='f_m',ylab='|Y_m|^2',col='red')
plot(fpad[(6000*3.7/12):(6000*4.3/12)+1],
     padpgram[(6000*3.7/12):(6000*4.3/12)+1],
     main='Filtered periodogram (fourth frequency interval)',
     xlab='f_m',ylab='|Y_m|^2',col='darkgreen')