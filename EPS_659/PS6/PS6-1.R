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
# Linear model
lm_co2 <- lm(co2~time)
s <- summary(lm_co2)
print(s)

# plot residuals vs. time
par(mar = c(5, 5, 5, 5))
plot(time,residuals(lm_co2),main="Residuals for Linear Model (ppm)",
     col="deeppink2", xlab="Time [years]", 
     ylab="Data Minus Model [ppm]",
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5) 

# plot QQ plot
qqnorm(residuals(lm_co2), main="QQ Plot for Linear Model",
     col="gold2", xlab="Model Values", ylab="Actual Values",
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)  
qqline(residuals(lm_co2), col="mediumseagreen", lty=2)

# get slope and y-intercept for linear model
slope <- s$coefficients[2,1] # concentration per month
cat(paste("Increase of CO2 concentration per month:", slope / 12))
cat(paste("Increase of CO2 concentration per year:",
          slope))
cat(paste("Increase of CO2 concentration per decade:",
          slope * 10))

# part (b)
print(paste("Linear model start-point:", predict(lm_co2)[1]))
cat(paste("Model minus actual at start-point:",
            predict(lm_co2)[1] - co2[1]))
cat(paste("Linear model endpoint:", 
            predict(lm_co2)[length(time)]))
cat(paste("Model minus actual at endpoint:", 
            predict(lm_co2)[length(time)] - co2[length(time)]))

# part (c)
y <- slope * time + s$coefficients[1,1]
std <- sd(y-co2)
cat(paste("St. dev. of the data-misfit residual (linear):",
            std))
MSM <- sum((y-mean(y))^2)
F <- MSM / (std ^ 2)
print(paste("F variance ratio:", F)) # F-ratio for linear model
1-pf(F, 1, length(co2)-1) # compute p-value 

# part (d)
# plot linear model against raw data
plot(time,co2,"l",main="Mauna Loa Carbon Dioxide (ppm),
     Linear Model",
     col="cyan3", xlab="Time [years]", 
     ylab="CO2 Per Unit Air [ppm]",
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
abline(lm_co2, col="maroon3")

# part (e)
# Quadratic model
time2 <- time^2
lm2_co2 <- lm(co2~time+time2)
s2 <- summary(lm2_co2)
print(s2)
pred2_co2 <- predict(lm2_co2)
print(paste("Quadratic model start-point:", pred2_co2[1]))
cat(paste("Model minus actual at start-point:",
          pred2_co2[1]-co2[1]))
cat(paste("Quadratic model endpoint:", 
          pred2_co2[length(time)]))
cat(paste("Model minus actual at endpoint:", 
          pred2_co2[length(time)]-co2[length(time)]))

# part (f)
cat(paste("Yearly increase in CO2 concentration per year", 
          "at the start of the time series:", 
          pred2_co2[13]-pred2_co2[1]))
cat(paste("Yearly increase in CO2 concentration per year", 
          "at the end of the time series:", 
          pred2_co2[762]-pred2_co2[750]))

# part (g)
y2 <- s2$coefficients[3,1] * time2 + 
  s2$coefficients[2,1] * time + s2$coefficients[1,1]
std2 <- sd(y2-co2)
cat(paste("St. dev. of the data-misfit residual (quadratic):", 
            std2))
MSM2 <- sum((y2-mean(y2))^2) 
# F-ratio for quadratic model
F2 <- (MSM2 / 2) / (sum((y2-co2)^2) / (length(time)-2))
print(paste("F variance ratio (quadratic case):", F2))
# p-value for quadratic model
1-pf(F2, 2, length(co2)-2)

# part (h)
# plot quadratic model against raw data
plot(time,co2,"l",main="Mauna Loa Carbon Dioxide (ppm),
     Quadratic Model",
     col="cyan3", xlab="Time [years]", 
     ylab="CO2 Per Unit Air [ppm]",
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
lines(time,y2, col="blue")

# plot residuals for quadratic model
plot(time,residuals(lm2_co2),main="Residuals for Quadratic Model (ppm)",
     col="magenta1", xlab="Time [years]", 
     ylab="Data Minus Model [ppm]",
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5) 

# run F-test to compare quadratic and linear models
Fcompare <- ((sum((y-co2)^2)-sum((y2-co2)^2)) / (2-1)) /
  (sum((y2-co2)^2) / (length(co2)-2))
print(paste("F variance ratio (comparing quadratic significance",
            "with linear case):", Fcompare))
# p-value for cubic model
1-pf(Fcompare, 2-1, length(co2)-2)

# part (i)
# Cubic model
time3 <- time^3
lm3_co2 <- lm(co2~time+time2+time3)
s3 <- summary(lm3_co2)
print(s3)
y3 <- s3$coefficients[4,1] * time3 + 
  s3$coefficients[3,1] * time2 + 
  s3$coefficients[2,1] * time + s3$coefficients[1,1]
std3 <- sd(y3-co2)
cat(paste("St. dev. of the data-misfit residual (cubic):",
          std3))

# run F-test to compare cubic and quadratic models
Fcompare <- ((sum((y2-co2)^2)-sum((y3-co2)^2)) / (3-2)) /
  (sum((y3-co2)^2) / (length(co2)-3))
print(paste("F variance ratio (comparing cubic significance",
"with quadratic case):", Fcompare))
# p-value for cubic model
1-pf(Fcompare, 3-2, length(co2)-3)

# plot cubic fit
plot(time,co2,"l",main="Mauna Loa Carbon Dioxide (ppm),
     Cubic Model",
     col="cyan3", xlab="Time [years]", 
     ylab="CO2 Per Unit Air [ppm]",
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
lines(time,y3, col="sienna1")

# plot cubic model residuals
plot(time,residuals(lm3_co2),main="Residuals for Cubic Model (ppm)",
     col="mediumaquamarine", xlab="Time [years]", 
     ylab="Data Minus Model [ppm]",
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5) 