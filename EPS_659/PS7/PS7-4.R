library(boot)
file <- paste0("C:/Users/jonas/Documents",
               "/EPS 659/PS7/co2_maunaloa.csv")
MaunaLoa <- read.csv(file)
year <- MaunaLoa$YEAR
month <- MaunaLoa$MONTH
time <- MaunaLoa$TIME
co2 <- MaunaLoa$CO2

model <- function(formula,data,indices){
  d <- data[indices,]
  fit <- lm(formula,data=d)
  return(coef(fit))
}

pi2 <- 2*pi
ann_s <- sin(pi2*time)
ann_c <- cos(pi2*time)
pi4 <- 4*pi
ann2_s <- sin(pi4*time)
ann2_c <- cos(pi4*time)
data_representer <- (time-1990.0)/25.0
data_representer2 <- data_representer^2
constant <- rep(1.0,length(data_representer))
gmatrix <- cbind(constant,data_representer,data_representer2,
                  ann_c,ann_s,ann2_c,ann2_s)
co2matrix <- data.frame(data=co2,g1=constant,g2=data_representer,
                        g3=data_representer2,g4=ann_c,
                        g5=ann_s,g6=ann2_c,g7=ann2_s)
results <- boot(data=co2matrix, statistic=model, R=1000,
                formula=data~g2+g3+g4+g5+g6+g7)
results
plot(results) # plot bootstrap results
 
# plot fit against data
plot(time,co2,"l",main="Mauna Loa Carbon Dioxide (ppm), 
     Bootstrapped Model",col="lightsalmon2", 
     xlab="Time [years]", 
     ylab="CO2 Per Unit Air [ppm]", 
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
lines(time,gmatrix %*% results$t0,col="darkolivegreen4")
grid()

# return coefficients and regression standard deviations
print("The bootstrap regression coefficients are")
print(results$t0)
print("The bootstrap regression standard deviations are")
print(apply(results$t, 2, sd))