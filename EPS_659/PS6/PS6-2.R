file <- paste0("C:/Users/jonas/Documents",
               "/EPS 659/PS6/co2_maunaloa.csv")
MaunaLoa <- read.csv(file)
names(MaunaLoa)
year <- MaunaLoa$YEAR
month <- MaunaLoa$MONTH
time <- MaunaLoa$TIME
rtime <- time-1990 
co2 <- MaunaLoa$CO2

# part (a)
# annual model (quadratic+annual cycle)
time2 <- time^2
pi2 <- 2*pi
ann_s <- sin(pi2*time)
ann_c <- cos(pi2*time)
lm2a_co2 <- lm(co2 ~ time + time2 + ann_s + ann_c)
s <- summary(lm2a_co2)
print(s)
pred <- coef(lm2a_co2)[[1]]+coef(lm2a_co2)[[2]]*time+
  coef(lm2a_co2)[[3]]*time2+
  coef(lm2a_co2)[[4]]*ann_s+coef(lm2a_co2)[[5]]*ann_c
# retrieve coefficients for annual cycle components
coef_s <- coef(lm2a_co2)[[4]]
coef_c <- coef(lm2a_co2)[[5]]
# retrieve errors for said coefficients
coef_se <- s$coefficients[[4,2]]
coef_ce <- s$coefficients[[5,2]]
# compute amplitude of annual cycle 
amp <- norm(c(coef_s, coef_c), type="2")
# compute error in amplitude of said cycle
ampe <- norm(c(coef_se, coef_ce), type="2")
cat(paste("The amplitude of the annual cycle is ",
          amp, "+/-", ampe))

# part (b)
# compute phase
phi <- atan2(coef_s,coef_c)
phi <- phi / pi2
cat(paste("The peak of the annual cycle is at", phi,
          "in terms of fraction of the year"))
cat(paste("The peak of the annual cycle is at", phi*12,
          "in terms of months"))


# quadratic fit
lm2_co2 <- lm(co2~time+time2)
s2 <- summary(lm2_co2)
pred2_co2 <- s2$coefficients[3,1] * time2 + 
  s2$coefficients[2,1] * time + s2$coefficients[1,1]

# best-fit cycle
cycle <- coef(lm2a_co2)[[4]]* sin(pi2*(1:100)*(1/100))+
  coef(lm2a_co2)[[5]]*cos(pi2*(1:100)*(1/100))
# cyclical temperature data arranged into months
arrangedinm <- matrix(co2[1:(12*floor(length(co2)/12))]-
                        pred2_co2[1:(12*floor(length(co2)/12))],
                      ncol=12, byrow=TRUE)
# annual cycle averaged over all years in data
cycledata <- colMeans(arrangedinm)
# plot annual cycle vs. average raw annual cycle
par(mar = c(5, 5, 5, 5))
plot((1:100)*(12/100),cycle,"l",main="Mauna Loa Averaged Annual 
Carbon Dioxide Cycle (ppm)",
     ylim=range(cycledata),
     col="dodgerblue2", xlab="Time [months]", 
     ylab="CO2 Per Unit Air [ppm]",
     cex.lab=1.5, cex.axis=1.0,
     cex.main=1.5, cex.sub=1.5)
points(3:12-0.5,cycledata[1:10], col="gold")
points(1:3-0.5,c(cycledata[11:12],cycledata[1]), 
      col="gold", lty=3)

# part (c)
# compute and plot residuals for quadratic vs. annual-cycle models
trendres <- co2 - pred2_co2
resres <- co2 - pred
plot(time,trendres,"l",main="Quadratic vs. Annual-Cycle 
     Model Residuals (ppm)",
     col="aquamarine3", xlab="Time [years]", 
     ylab="CO2 Per Unit Air [ppm]",
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
lines(time,resres, col="darkorange1")

cat(paste("The quadratic model residual st. dev. is",
            sd(trendres)))
cat(paste("The quadratic+annual-cycle model residual",
          "st. dev. is", sd(resres)))

# part (d)
# compute quadratic+annual cycle+semi-annual cycle model
pi4 <- 4*pi
ann2_s <- sin(pi4*time)
ann2_c <- cos(pi4*time)
lm2a2_co2 <- lm(co2 ~ time + time2 + ann_s + 
                  ann_c + ann2_s + ann2_c)
sa2 <- summary(lm2a2_co2)
print(sa2)
# retrieve coefficients for semi-annual cycle components
coef_s2 <- coef(lm2a2_co2)[[6]]
coef_c2 <- coef(lm2a2_co2)[[7]]
# retrieve errors for said coefficients
coef_s2e <- sa2$coefficients[[6,2]]
coef_c2e <- sa2$coefficients[[7,2]]
# compute amplitude of semi-annual cycle 
amp2 <- norm(c(coef_s2, coef_c2), type="2")
# compute error in said cycle
amp2e <- norm(c(coef_s2e, coef_c2e), type="2")
cat(paste("The amplitude of the semi-annual cycle is ",
          amp2, "+/-", amp2e))

# best-fit annual+semi-annual cycle
cycle2 <- coef(lm2a2_co2)[[4]]* sin(pi2*(1:100)*(1/100))+
  coef(lm2a2_co2)[[5]]*cos(pi2*(1:100)*(1/100))+
  coef(lm2a2_co2)[[6]]* sin(pi4*(1:100)*(1/100))+
  coef(lm2a2_co2)[[7]]*cos(pi4*(1:100)*(1/100))
# plot annual+semi-annual cycle vs. annual only cycle vs. data
par(mar = c(5, 5, 5, 5))
plot((1:100)*(12/100),cycle2,"l",main="Mauna Loa Averaged Annual 
Carbon Dioxide Cycle (ppm)",
     ylim=range(cycledata),
     col="darkolivegreen3", xlab="Time [months]", 
     ylab="CO2 Per Unit Air [ppm]",
     cex.lab=1.5, cex.axis=1.0,
     cex.main=1.5, cex.sub=1.5)
lines((1:100)*(12/100),cycle,col="dodgerblue2")
points(3:12-0.5,cycledata[1:10], col="gold")
points(1:3-0.5,c(cycledata[11:12],cycledata[1]), 
       col="gold", lty=3)

# part (e)
# full residual for quadratic+annual+semi-annual model
newres <- co2 - (coef(lm2a2_co2)[[1]]+coef(lm2a2_co2)[[2]]*time+
                   coef(lm2a2_co2)[[3]]*time2+
                   coef(lm2a2_co2)[[4]]*ann_s+
                   coef(lm2a2_co2)[[5]]*ann_c+
                   coef(lm2a2_co2)[[6]]*ann2_s+
                   coef(lm2a2_co2)[[7]]*ann2_c)
# plot residuals for said model
plot(time,newres,"l",main="Final Model Residuals (ppm)",
     col="darkorchid3", xlab="Time [years]", 
     ylab="CO2 Per Unit Air [ppm]",
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)

cat(paste("The quadratic+annual+semi-annual-cycle model residual",
          "st. dev. is", sd(newres)))

# plot of final data model
yfull <- coef(lm2a2_co2)[[1]]+coef(lm2a2_co2)[[2]]*time+
  coef(lm2a2_co2)[[3]]*time2+
  coef(lm2a2_co2)[[4]]*ann_s+
  coef(lm2a2_co2)[[5]]*ann_c+
  coef(lm2a2_co2)[[6]]*ann2_s+
  coef(lm2a2_co2)[[7]]*ann2_c
plot(time,co2,"l",main="Mauna Loa Carbon Dioxide (ppm),
     Full Model",
     col="cyan3", xlab="Time [years]", 
     ylab="CO2 Per Unit Air [ppm]",
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
lines(time,yfull, col="sienna1")