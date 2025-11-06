file <- paste0("C:/Users/jonas/Documents",
               "/EPS 659/PS7/co2_maunaloa.csv")
MaunaLoa <- read.csv(file)
# names(MaunaLoa)
year <- MaunaLoa$YEAR
month <- MaunaLoa$MONTH
time <- MaunaLoa$TIME
co2 <- MaunaLoa$CO2
# plot(time,co2,"l")
# abline(lm(co2~time))
# title(main="Mauna Loa Carbon Dioxide (ppm)")

time2 <- time^2
pi2 <- 2*pi
ann_s <- sin(pi2*time)
ann_c <- cos(pi2*time)
pi4 <- 4*pi
ann2_s <- sin(pi4*time)
ann2_c <- cos(pi4*time)
# pre-define arrays for regression parameters 
# of seven data segments
ann1 <- 1:7 # for part (a)
ann1_err <- 1:7
ann_phase <- 1:7 # for part (b)
ann_pherr <- 1:7
ann2 <- 1:7 # for part (c)
ann2_err <- 1:7
yr <- 1:7
# initialize segment counter
k <- 0
# the data segments are defined by their start-points
for (j in (1:7)*120-119){
  k <- k+1
  timed <- time[j:(j+119)]
  time2d <- timed^2
  co2d <- co2[j:(j+119)]
  ann_sd <- ann_s[j:(j+119)]
  ann_cd <- ann_c[j:(j+119)]
  ann2_sd <- ann2_s[j:(j+119)]
  ann2_cd <- ann2_c[j:(j+119)]
  lmlmlm <- lm(co2d ~ timed + time2d + ann_sd +
                 ann_cd + ann2_sd + ann2_cd)
  anns <- coef(summary(lmlmlm))["ann_sd","Estimate"]
  annc <- coef(summary(lmlmlm))["ann_cd","Estimate"]
  ann_amp <- sqrt(anns^2 + annc^2)
  ann1[k] <- ann_amp 
  ann_phase[k] <- 180*atan2(anns,annc)/pi
  yr[k] <- year[j+60]
  anns <- coef(summary(lmlmlm))["ann_sd","Std. Error"]
  annc <- coef(summary(lmlmlm))["ann_cd","Std. Error"]
  ann1_err[k] <- sqrt(anns^2+annc^2)
  # the phase uncertainty in radians is the relative 
  # uncertainty in coefficient amplitude
  # multiply this phase uncertainty by 180/pi to express 
  # uncertainty in degrees
  ann_pherr[k] <- (180/pi)*sqrt(anns^2+annc^2)/ann_amp
  ann2s <- coef(summary(lmlmlm))["ann2_sd","Estimate"]
  ann2c <- coef(summary(lmlmlm))["ann2_cd","Estimate"]
  ann2[k] <- sqrt(ann2s^2 + ann2c^2)
  ann2s <- coef(summary(lmlmlm))["ann2_sd","Std. Error"]
  ann2c <- coef(summary(lmlmlm))["ann2_cd","Std. Error"]
  ann2_err[k] <- sqrt(ann2s^2+ann2c^2)
}

# part (a)
# plot annual-cycle amplitudes
par(mar = c(5, 5, 5, 5))
plot(yr,ann1,"l",main="Mauna Loa Carbon Dioxide (ppm),
     Annual-Cycle Amplitudes",
     col="mediumorchid3", xlab="Decade [years]", 
     ylab="CO2 Per Unit Air [ppm]", 
     ylim=c(min(ann1-ann1_err),max(ann1+ann1_err)),
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
arrows(x0=yr,y0=ann1+ann1_err,x1=yr,y1=ann1-ann1_err,
       code=3,angle=90,length=0.2)
grid()
print("The annual-cycle amplitudes are")
cat(paste(ann1,"+/-",ann1_err))

# compute F-test to see if amplitude change is significant
# explained variance
expvar <- 120 * sum((ann1-mean(ann1))^2) / (7-1) 
# unexplained variance
unexpvar <- sum(ann1_err^2)*(120-(1+6))/(120*7-7) 
Fratio <- expvar / unexpvar
# p-value for changing amplitude vs. constant
1-pf(Fratio, 7-1, 120*7-7)

# part (b)
# plot phase
plot(yr,ann_phase,"l",main="Mauna Loa Carbon Dioxide (ppm),
     Annual-Cycle Phases",
     col="lightsalmon2", xlab="Decade [years]", 
     ylab="Degrees", 
     ylim=c(min(ann_phase-ann_pherr),max(ann_phase+ann_pherr)),
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
arrows(x0=yr,y0=ann_phase+ann_pherr,x1=yr,y1=ann_phase-ann_pherr,
       code=3,angle=90,length=0.2)
grid()
print("The annual-cycle phases are")
cat(paste(ann_phase,"+/-",ann_pherr))

# compute F-test to see if phase change is significant
# explained variance
expvar <- 120 * sum((ann_phase-mean(ann_phase))^2) / (7-1) 
# unexplained variance
unexpvar <- sum(ann_pherr^2)*(120-(1+6))/(120*7-7) 
Fratio <- expvar / unexpvar
# p-value for changing amplitude vs. constant
1-pf(Fratio, 7-1, 120*7-7)

# part (c)
# plot semi-annual-cycle amplitudes
plot(yr,ann2,"l",main="Mauna Loa Carbon Dioxide (ppm),
     Semi-Annual-Cycle Amplitudes",
     col="olivedrab3", xlab="Decade [years]", 
     ylab="CO2 Per Unit Air [ppm]", 
     ylim=c(min(ann2-ann2_err),max(ann2+ann2_err)),
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
arrows(x0=yr,y0=ann2+ann2_err,x1=yr,y1=ann2-ann2_err,
       code=3,angle=90,length=0.2)
grid()
print("The semi-annual-cycle amplitudes are")
cat(paste(ann2,"+/-",ann2_err))

# compute F-test to see if amplitude change is significant
# explained variance
expvar <- 120 * sum((ann2-mean(ann2))^2) / (7-1) 
# unexplained variance
unexpvar <- sum(ann2_err^2)*(120-(1+6))/(120*7-7) 
Fratio <- expvar / unexpvar
# p-value for changing amplitude vs. constant
1-pf(Fratio, 7-1, 120*7-7)