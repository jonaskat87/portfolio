file <- paste0("C:/Users/jonas/Documents",
               "/EPS 659/PS7/co2_maunaloa.csv")
MaunaLoa <- read.csv(file)
year <- MaunaLoa$YEAR
month <- MaunaLoa$MONTH
time <- MaunaLoa$TIME
co2 <- MaunaLoa$CO2

# part (a)
em <- lm(log(co2)~time) # coefficients are log(B), C
Btilde <- coef(summary(em))["(Intercept)","Estimate"]
C <- coef(summary(em))["time","Estimate"]
summary(em)

# plot exponential model vs. data (semi-log plot)
par(mar = c(5, 5, 5, 5))
plot(time,Btilde + C*time,"l",
main="Mauna Loa Carbon Dioxide (ppm),
     Exponential Model", ylim=range(log(co2)),
     col="lightseagreen", xlab="Time [years]", 
     ylab="Logarithm of CO2 Per Unit Air [log(ppm)]", 
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
lines(time,log(co2), col="lightsalmon2")
grid()
plot(time,log(co2)-Btilde-C*time,"l",
main="Mauna Loa Carbon Dioxide (ppm),
     Exponential Model Residuals",
     col="palevioletred3", xlab="Time [years]", 
     ylab="Logarithm of CO2 Per Unit Air [log(ppm)]", 
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
grid()

# part (b)
# plot exponential model vs. data (physical plot)
plot(time,exp(Btilde+C*time),"l",
     main="Mauna Loa Carbon Dioxide (ppm),
     Exponential Model", ylim=range(co2),
     col="lightseagreen", xlab="Time [years]", 
     ylab="CO2 Per Unit Air [ppm]", 
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
lines(time,co2, col="lightsalmon2")
grid()
plot(time,co2-exp(Btilde+C*time),"l",
     main="Mauna Loa Carbon Dioxide (ppm),
     Exponential Model Residuals",
     col="palevioletred3", xlab="Time [years]", 
     ylab="CO2 Per Unit Air [ppm]", 
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
grid()
cat(paste("The endpoints of the model are", 
          exp(Btilde+C*time[1]), "and", 
          tail(exp(Btilde+C*time), n=1), "which is 
          associated with residuals", co2[1]-
            exp(Btilde+C*time[1]), "and",
          tail(co2-exp(Btilde+C*time), n=1), 
          "respectively."))

# part (c)
# nonlinear fit of log(co2)
# subtract 1788 from time before regression, 
# to use pre-industrial initial condition the nls function 
# in R asks for initial values for the regression parameters
# which are identified by the start parameter
time0 <- time-1788
# use parameters above with A=residual from part (b) at 
# start of data for an initial guess
astart <- co2[1]-exp(Btilde+C*time[1])
bstart <- exp(Btilde)
cstart <- C
nlm_co2 <- nls(co2~a+b*exp(c*time0), alg = "port",
               control=nls.control(maxiter = 500),
               start=list(a=astart,b=bstart,c=cstart))
predexp <- predict(nlm_co2)
summary(nlm_co2)
a <- coef(summary(nlm_co2))["a","Estimate"]
b <- coef(summary(nlm_co2))["b","Estimate"]
c <- coef(summary(nlm_co2))["c","Estimate"]
modelval <- a+b*exp(c*time0)

plot(time,co2,"l",main="Mauna Loa Carbon Dioxide (ppm), 
     Nonlinear Model",col="lightsalmon2", 
     xlab="Time [years]", 
     ylab="CO2 Per Unit Air [ppm]", 
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
lines(time,modelval,col="purple3")
grid()
plot(time,co2-modelval,"l",
     main="Mauna Loa Carbon Dioxide (ppm),
     Nonlinear Model Residuals",
     col="blue3", xlab="Time [years]", 
     ylab="CO2 Per Unit Air [ppm]", 
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
grid()

time2 <- time^2
lm2_co2 <- lm(co2~time+time2)
s2 <- summary(lm2_co2)
y2 <- s2$coefficients[3,1] * time2 + 
  s2$coefficients[2,1] * time + s2$coefficients[1,1]
y2var <- sum((y2-co2)*(y2-co2))
nlvar <- sum((modelval-co2)*(modelval-co2))
cat(paste("The residual variance of our nonlinear model is",
          nlvar))
cat(paste("The residual variance of the quadratic model is",
          y2var))
cat(paste("St. dev. of the data-misfit residual (nonlinear):", 
          sd(modelval-co2)))
cat(paste("St. dev. of the data-misfit residual (quadratic):", 
          sd(y2-co2)))
# run F-test to compare nonlinear and quadratic models
Fcompare <- -((nlvar - y2var) / (3-2)) /
  (y2var / (length(co2)-3))
cat(paste("F variance ratio (comparing nonlinear model 
            significance with linear case):", Fcompare))
# p-value for nonlinear vs. quadratic model
1-pf(Fcompare, 3-2, length(co2)-3)

# part (d)
# standard deviation for the first term
aerr <- coef(summary(nlm_co2))["a","Std. Error"]
cat(paste("Our model's prediction for the preindustrial level
          of CO2 is", a, "+/-", aerr))
cat(paste("Our computed value for A is between", 
          abs((a-260)/aerr), "and", abs((a-270)/aerr), 
          "standard deviations from the official estimates"))

# part (e)
berr <- coef(summary(nlm_co2))["b","Std. Error"]
cat(paste("Our model's prediction for the year 1788 is", a+b,
            "+/-", sqrt(aerr^2 + berr^2)))