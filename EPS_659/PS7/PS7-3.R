file <- paste0("C:/Users/jonas/Documents",
               "/EPS 659/PS7/co2_maunaloa.csv")
MaunaLoa <- read.csv(file)
year <- MaunaLoa$YEAR
month <- MaunaLoa$MONTH
time <- MaunaLoa$TIME
co2 <- MaunaLoa$CO2

# part (a)
pi2 <- 2*pi
ann_s <- sin(pi2*time)
ann_c <- cos(pi2*time)
pi4 <- 4*pi
ann2_s <- sin(pi4*time)
ann2_c <- cos(pi4*time)
data_representer <- (time-1990.0)/25.0
data_representer2 <- data_representer^2
constant <- rep(1.0,length(data_representer))
gmatrixa <- cbind(constant,data_representer,data_representer2,
                 ann_c,ann_s,ann2_c,ann2_s)
svddataa <- svd(gmatrixa)
k <- 1 : 7

# plot part (a) singular values against k
plot(k,svddataa$d,pch=16,
     main="Singular Values of the Data Representation Matrix, G
     (Part A)",
     col="navyblue", xlab="k", 
     ylab="Singular Value [-]", 
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
cat(paste("The condition number of gmatrix (part a) is", 
          svddataa$d[1] / svddataa$d[7]))
grid()

# plot non-periodic data representers for part (a)
plot(time,constant,type='l',
     ylim=range(c(data_representer,data_representer2)),
     main="Non-Periodic Data Representers for Part A",
     col="red3", xlab="Time [years]", 
     ylab="Singular Value [-]", 
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
lines(time,data_representer, col="mediumblue")
lines(time,data_representer2, col="springgreen3")
grid()

# part (b)
print("The principal component of the model space with the
      largest singular value is")
print(svddataa$v[,1])
print("The principal component of the model space with the
      smallest singular value is")
print(svddataa$v[,7])
print("The dot product between these components is")
print(svddataa$v[,1]%*%svddataa$v[,7])

# part (c)
data_representer <- time/2000.0
data_representer2 <- data_representer^2
gmatrixc <- cbind(constant,data_representer,data_representer2,
                 ann_c,ann_s,ann2_c,ann2_s)
svddatac <- svd(gmatrixc)
k <- 1 : 7

# plot part (c) singular values against k
plot(k,svddatac$d,pch=16,
     main="Singular Values of the Data Representation Matrix, G
     (Part C)",
     col="orangered3", xlab="k", 
     ylab="Singular Value [-]", 
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
cat(paste("The condition number of gmatrix (part c) is", 
          svddatac$d[1] / svddatac$d[7]))
grid()

# plot non-periodic data representers for part (c)
plot(time,constant,type='l',
     ylim=range(c(data_representer,data_representer2)),
     main="Non-Periodic Data Representers for Part C",
     col="red3", xlab="Time [years]", 
     ylab="Singular Value [-]", 
     cex.lab=1.5, cex.axis=1.0, cex.main=1.5, cex.sub=1.5)
lines(time,data_representer, col="mediumblue")
lines(time,data_representer2, col="springgreen3")
grid()