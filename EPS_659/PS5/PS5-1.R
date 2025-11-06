fileglobal <- paste0("C:/Users/jonas/Documents",
              "/EPS 659/PS5/Hadcrut_GlobalAverage.csv")
had_ns_avg <- read.csv(fileglobal)
temps <- had_ns_avg$DTEMPC
years <- had_ns_avg$YEAR
months <- had_ns_avg$MONTH

time <- years + (months-0.5)/12.0
fileSH <- paste0("C:/Users/jonas/Documents",
                    "/EPS 659/PS5/Hadcrut_SHAverage.csv")
fileNH <- paste0("C:/Users/jonas/Documents",
                "/EPS 659/PS5/Hadcrut_NHAverage.csv")
had_sh_avg <- read.csv(fileSH)
had_nh_avg <- read.csv(fileNH)
temps_sh <- had_sh_avg$DTEMPC
temps_nh <- had_nh_avg$DTEMPC

# part (a)
# time is indexed in terms of months. 
# Hence, in terms of years, t=((i-0.5)/12)+1850
# The indices for the first decade of the 20th-century
# start when 1900=((i-0.5)/12)+1850, such that i=601
# Similarly, the indices for the last decade of the 20th-century
# start when 1990=((i-0.5)/12)+1850, such that i=1681
firstmean <- mean(temps[1:(12*10)+600]) # sample mean of 1900s temp
firstsd <- sd(temps[1:(12*10)+600]) # sample st. dev. of 1900s temp
lastmean <- mean(temps[1:(12*10)+1680]) # sample mean of 1990s temp
lastsd <- sd(temps[1:(12*10)+1680]) # sample st. dev. of 1990s temp
print(paste("sample mean for the first decade:", firstmean))
print(paste("st. dev. for the first decade:", firstsd))
print(paste("sample mean for the last decade:", lastmean))
print(paste("st. dev. for the last decade:", lastsd))
tfirst <- sqrt(12*10-1)*firstmean / firstsd # t for first decade
tlast <-  sqrt(12*10-1)*lastmean / lastsd # t for last decade
print(paste("t for the first decade:", tfirst))
print(paste("Probability of type 1 error:", pt(-abs(tfirst), 12*10-1)*2))
print(paste("t for the last decade:", tlast))
print(paste("Probability of type 1 error:", pt(-abs(tlast), 12*10-1)*2))

# part (b)
# two-way t-test
twot <- sqrt(12*10-1)*(firstmean-lastmean) / sqrt(firstsd^2+lastsd^2)
print(paste("t between first and last decades:", twot))
print(paste("Probability of type 1 error:", pt(-abs(twot), 2*(12*10-1))*2))

# part (c)
# index vector for indices at beginning of each decade
indices <- seq(0,length(temps)-12*10,by=12*10)
tval <- rep(0,length(indices)-1) # vector to store t-values in
# loop through adjacent decades and compute t-values
for (i in 1:(length(indices)-1)) {
  firstmean <- mean(temps[1:(12*10)+indices[i]])
  firstsd <- sd(temps[1:(12*10)+indices[i]]) 
  lastmean <- mean(temps[1:(12*10)+indices[i+1]]) 
  lastsd <- sd(temps[1:(12*10)++indices[i+1]]) 
  tfirst <- sqrt(12*10-1)*firstmean / firstsd 
  tlast <-  sqrt(12*10-1)*lastmean / lastsd 
  tval[i] <- sqrt(12*10-1)*(firstmean-lastmean) / 
    sqrt(firstsd^2+lastsd^2)
}
# plot probability of type 1 errors vs. decade pair
par(mar = c(5, 5, 5, 5))
plot(time[indices]-10,pt(-abs(tval),2*(12*10-1))*2,log="y",
     main="Type 1 Error Probabilities in the UK Hadley Center
 Global Average Data Series from 1850-2010",
     xlab="Earlier Decade in Decadal Pair [years]",
     ylab="Probability of Type 1 Error", col="purple", pch=19, 
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
grid()
print(paste("Probability of type 1 error:", pt(-abs(tval),2*(12*10-1))*2))

# test increases vs. decreases at 99% confidence level
# test data for if we reject null hypothesis (1=yes, 0=no)
# alternative: first decade is greater than second
morepval <- rep(0,length(indices)-1) # p-values for decrease
# alternative: first decade is less than second
lesspval <- rep(0,length(indices)-1) # p-values for increase
for (i in 1:(length(indices)-1)) {
  ttestmore <- t.test(temps[1:(12*10)+indices[i]],
                 temps[1:(12*10)+indices[i+1]],
                 alternative="greater", conf.level = 0.99)
  ttestless <- t.test(temps[1:(12*10)+indices[i]],
                        temps[1:(12*10)+indices[i+1]],
                        alternative="less", conf.level = 0.99)
  morepval[i] <- ttestmore$p.value
  lesspval[i] <- ttestless$p.value
}
cat("P-values for null hypothesis of no significant increase
      in temperature globally:")
print(lesspval)
plot(time[indices]-10,lesspval,log="y",
     main="P-values for Rejecting Null Hypothesis 
     w.r.t Increase (Globally)",
     xlab="Earlier Decade in Decadal Pair [years]",
     ylab="P-Value", col="purple", pch=19, 
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
grid()
cat("P-values for null hypothesis of no significant decrease
      in temperature globally:")
print(morepval)
plot(time[indices]-10,morepval,log="y",
     main="P-values for Rejecting Null Hypothesis 
     w.r.t Decrease (Globally)",
     xlab="Earlier Decade in Decadal Pair [years]",
     ylab="P-Value", col="purple", pch=19, 
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
grid()

# part (d)
tvalN <- rep(0,length(indices)-1) 
for (i in 1:(length(indices)-1)) {
  firstmean <- mean(temps_nh[1:(12*10)+indices[i]])
  firstsd <- sd(temps_nh[1:(12*10)+indices[i]]) 
  lastmean <- mean(temps_nh[1:(12*10)+indices[i+1]]) 
  lastsd <- sd(temps_nh[1:(12*10)+indices[i+1]]) 
  tfirst <- sqrt(12*10-1)*firstmean / firstsd 
  tlast <-  sqrt(12*10-1)*lastmean / lastsd 
  tvalN[i] <- sqrt(12*10-1)*(firstmean-lastmean) / 
    sqrt(firstsd^2+lastsd^2)
}

par(mar = c(5, 5, 5, 5))
plot(time[indices]-10,pt(-abs(tvalN),2*(12*10-1))*2,log="y",
     main="Type 1 Error Probabilities in the UK Hadley Center
 Northern Hemisphere Average Data Series from 1850-2010",
     xlab="Earlier Decade in Decadal Pair [years]",
     ylab="Probability of Type 1 Error", col="red", pch=19, 
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
grid()
print(paste("Probability of type 1 error:", pt(-abs(tvalN),2*(12*10-1))*2))

morepvalN <- rep(0,length(indices)-1) 
lesspvalN <- rep(0,length(indices)-1) 
for (i in 1:(length(indices)-1)) {
  ttestmore <- t.test(temps_nh[1:(12*10)+indices[i]],
                     temps_nh[1:(12*10)+indices[i+1]],
                     alternative="greater", conf.level = 0.99)
  ttestless <- t.test(temps_nh[1:(12*10)+indices[i]],
                     temps_nh[1:(12*10)+indices[i+1]],
                     alternative="less", conf.level = 0.99)
  morepvalN[i] <- ttestmore$p.value
  lesspvalN[i] <- ttestless$p.value
}
cat("P-values for null hypothesis of no significant increase
      in temperature in the Northern Hemisphere:")
print(lesspvalN)
plot(time[indices]-10,lesspvalN,log="y",
     main="P-values for Rejecting Null Hypothesis 
     w.r.t Increase (Northern Hemisphere)",
     xlab="Earlier Decade in Decadal Pair [years]",
     ylab="P-Value", col="red", pch=19, 
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
grid()
cat("P-values for null hypothesis of no significant decrease
      in temperature in the Northern Hemisphere:")
print(morepvalN)
plot(time[indices]-10,morepvalN,log="y",
     main="P-values for Rejecting Null Hypothesis 
     w.r.t Decrease (Northern Hemisphere)",
     xlab="Earlier Decade in Decadal Pair [years]",
     ylab="P-Value", col="red", pch=19, 
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
grid()

# part (e)
tvalS <- rep(0,length(indices)-1) 
pvalS <- rep(0,length(indices)-1)
for (i in 1:(length(indices)-1)) {
  firstmean <- mean(temps_sh[1:(12*10)+indices[i]])
  firstsd <- sd(temps_sh[1:(12*10)+indices[i]]) 
  lastmean <- mean(temps_sh[1:(12*10)+indices[i+1]]) 
  lastsd <- sd(temps_sh[1:(12*10)+indices[i+1]]) 
  tfirst <- sqrt(12*10-1)*firstmean / firstsd 
  tlast <-  sqrt(12*10-1)*lastmean / lastsd 
  tvalS[i] <- sqrt(12*10-1)*(firstmean-lastmean) / 
    sqrt(firstsd^2+lastsd^2)
}

par(mar = c(5, 5, 5, 5))
plot(time[indices]-10,pt(-abs(tvalS),2*(12*10-1))*2,log="y",
     main="Type 1 Error Probabilities in the UK Hadley Center
 Southern Hemisphere Average Data Series from 1850-2010",
     xlab="Earlier Decade in Decadal Pair [years]",
     ylab="Probability of Type 1 Error", col="blue", pch=19, 
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
grid()
print(paste("Probability of type 1 error:", pt(-abs(tvalS),2*(12*10-1))*2))

morepvalS <- rep(0,length(indices)-1) 
lesspvalS <- rep(0,length(indices)-1) 
for (i in 1:(length(indices)-1)) {
  ttestmore <- t.test(temps_sh[1:(12*10)+indices[i]],
                     temps_sh[1:(12*10)+indices[i+1]],
                     alternative="greater", conf.level = 0.99)
  ttestless <- t.test(temps_sh[1:(12*10)+indices[i]],
                     temps_sh[1:(12*10)+indices[i+1]],
                     alternative="less", conf.level = 0.99)
  morepvalS[i] <- ttestmore$p.value
  lesspvalS[i] <- ttestless$p.value
}
cat("P-values for null hypothesis of no significant increase
      in temperature in the Southern Hemisphere:")
print(lesspvalS)
plot(time[indices]-10,lesspvalS,log="y",
     main="P-values for Rejecting Null Hypothesis 
     w.r.t Increase (Southern Hemisphere)",
     xlab="Earlier Decade in Decadal Pair [years]",
     ylab="P-Value", col="blue", pch=19, 
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
grid()
cat("P-values for null hypothesis of no significant decrease
    in temperature in the Southern Hemisphere:")
print(morepvalS)
plot(time[indices]-10,morepvalS,log="y",
     main="P-values for Rejecting Null Hypothesis 
     w.r.t Decrease (Southern Hemisphere)",
     xlab="Earlier Decade in Decadal Pair [years]",
     ylab="P-Value", col="blue", pch=19, 
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
grid()