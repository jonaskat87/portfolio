file = paste0("C:/Users/jonas/Documents",
              "/EPS 659/PS4/GlobalEarthquakes6.0.csv")
gevents <- read.csv(file)
years <- gevents$YEAR
months <- gevents$MONTH
days <- gevents$DAY
mag <- gevents$MAG

# part (c)
sortedmag <- sort(mag)
ecdfeval <- ecdf(sortedmag)
eCDF <- ecdfeval(sortedmag) 
par(mar = c(7, 7, 7, 7))
plot(sortedmag, eCDF, type="l", col="brown", 
     main="Global Earthquakes with Richter Magnitudes M >= 6.0
     from 1980 to 2020 (Data from IRS)",
     sub=" (Data from IRS)",
     xlab="Richter Magnitude", 
     ylab="CDF, F(M) 
     (Proportion with Magnitude <= M)",
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
grid()

# part (d)
par(mar = c(5, 5, 5, 5))
plot(1-eCDF, sortedmag, type="l", col="brown", log="y",
     main="Global Earthquakes with Richter Magnitudes M >= 6.0
     from 1980 to 2020",
     sub=" (Data from IRS)",
     xlab="1-F(M) (Proportion with Magnitude >= M)",
     ylab="Richter Magnitude", 
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
grid()

# part (e)
# for counting the number of data values in a given year
count8 <- rep(0,length(1980:2020))
# Retrieve years corresponding to earthquakes with M>=8.0
years8 <- years[mag>=8] 
# we want a function that, given a year, will return the index
# in count corresponding to that year. Clearly, this function 
# can be linear, since the years chosen (T) follow a linear 
# trend and so do the indices that correspond to those years.
# The function is i=f(T)=T-1979.
# and now we count earthquakes for each year!
for (i in 1:length(years8)) {
  count8[years8[i]-1979] <- count8[years8[i]-1979] + 1
}
lambda8 <- mean(count8) # sample mean
# for counting the number of times we have a certain 
# number of earthquakes each year (starting with 0)
occur8 <- rep(0,max(count8)) 
# count number of occurrences for given counts/year
for (i in 1:length(occur8)) {
  occur8[i] = sum(count8==(i-1))
}
occur8 <- occur8 / sum(occur8) # normalize PDF

CDF8 <- ppois(-1:6, lambda8) # compute Poisson CDF
# plot CDF
plot(-1:6, CDF8, type="p", col="darkgoldenrod", 
     main=paste("Poisson Distribution CDF with lambda=", lambda8),
     xlab="K=k (Number of Events)",
     ylab="CDF", pch=19,
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
grid()

PDF8 <- dpois(-1:6, lambda8) # compute Poisson PDF
# plot Poisson PDF against empirical PDF
plot(-1:6, PDF8, type="p", col="darkgoldenrod", 
     main=paste("Poisson Distribution PDF with lambda=", lambda8, "
                vs. Empirical PDF"),
     ylim=range(c(PDF8, occur8)),
     xlab="K=k (Number of Events)",
     ylab="PDF", pch=19,
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
lines(-1:(length(occur8)+2), c(0,occur8,rep(0,3)), type="p",
      col="firebrick", pch=19)
grid()
# compute relative L1 error between PDFs
L1error8 = sum(abs(c(0,occur8,rep(0,3))-PDF8)) / sum(abs(PDF8))
print(paste("The relative L1 error for the M>=8.0 case is", L1error8))
# compute relative L2 error between PDFs
L2error8 = sum((c(0,occur8,rep(0,3))-PDF8)^2) / sum((PDF8)^2)
print(paste("The relative L2 error for the M>=8.0 case is", L2error8))

# part (f)
count7.5 <- rep(0,length(1980:2020))
years7.5 <- years[mag>=7.5] 
for (i in 1:length(years7.5)) {
  count7.5[years7.5[i]-1979] <- count7.5[years7.5[i]-1979] + 1
}
lambda7.5 <- mean(count7.5) 

occur7.5 <- rep(0,max(count7.5)) 
for (i in 1:length(occur7.5)) {
  occur7.5[i] = sum(count7.5==(i-1))
}
occur7.5 <- occur7.5 / sum(occur7.5) 
CDF7.5 <- ppois(-1:12, lambda7.5) 
plot(-1:12, CDF7.5, type="p", col="darkgoldenrod", 
     main=paste("Poisson Distribution CDF with lambda=", lambda7.5),
     xlab="K=k (Number of Events)",
     ylab="CDF", pch=19,
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
grid()

PDF7.5 <- dpois(-1:12, lambda7.5) 
plot(-1:12, PDF7.5, type="p", col="darkgoldenrod", 
     main=paste("Poisson Distribution PDF with lambda=", lambda7.5,
                "vs. Empirical PDF"),
     ylim=range(c(PDF7.5, occur7.5)),
     xlab="K=k (Number of Events)",
     ylab="PDF", pch=19,
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
lines(-1:(length(occur7.5)+2), c(0,occur7.5,rep(0,3)), type="p",
      col="firebrick", pch=19)
grid()

L1error7.5 = sum(abs(c(0,occur7.5,rep(0,3))-PDF7.5)) / sum(abs(PDF7.5))
print(paste("The relative L1 error for the M>=7.5 case is", L1error7.5))

L2error7.5 = sum((c(0,occur7.5,rep(0,3))-PDF7.5)^2) / sum((PDF7.5)^2)
print(paste("The relative L2 error for the M>=7.5 case is", L2error7.5))