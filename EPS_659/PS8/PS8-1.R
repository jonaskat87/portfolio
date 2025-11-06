# part (b)
dx <- 0.01 # grid size for plot
Gamma0 <- 6.67408e-11
# define data representer via a function
g <- function(x,xi,h,Dz) {2*Gamma0*h*Dz/(((x-xi)^2+h^2))
}
hq <- (1:3)*5*1000
Dz <- 1000
N <- ((80+60)/dx)+1
xq <- (((0:(N-1))*dx)-60)*1000

# plot data representers
len <- length(xq)
par(mar=c(5,5,5,5))
all <- c(g(0,xq,hq[1],Dz),g(10*1000,xq,hq[1],Dz),
         g(0,xq,hq[2],Dz),g(10*1000,xq,hq[2],Dz),
         g(0,xq,hq[3],Dz),g(10*1000,xq,hq[3],Dz))
plot(xq,all[1:len],type='l',col='red',xaxt="n",
     xlab='xi [km]', main='Sample data representers',
     ylim=range(all), ylab='g(x,xi) [m^(3)kg^(-1)s^(-2)]')
lines(xq,all[(1:len)+len],'l',col='gold',xaxt="n")
lines(xq,all[(1:len)+2*len],'l',col='orange',xaxt="n")
lines(xq,all[(1:len)+3*len],'l',col='green',xaxt="n")
lines(xq,all[(1:len)+4*len],'l',col='blue',xaxt="n")
lines(xq,all[(1:len)+5*len],'l',col='purple',xaxt="n")
axis(1, at=((0:7)*20-60)*1000, labels=(0:7)*20-60)

# part (c)
file <- paste0("C:/Users/jonas/Documents",
               "/EPS 659/PS8/gravity_pset8.csv")
gdata <- read.csv(file)
X <- gdata$X
n <- length(X)

hq <- c(5,10,20)*1000
const <- (pi*(2*Gamma0*Dz)^2)/(2*hq)
diff <- outer(X, X, "-") 
eigval <- matrix(0*(1:(3*n)), nrow = 3) # to store eigenvalues
cond <- 1:3 # to store condition numbers
for (i in 1:3)
{
  Gamma <- const[i] / ((diff/(2*hq[i]))^2+1)
  cond[i] <- norm(Gamma,'2')*norm(solve(Gamma),'2')
  ev <- eigen(Gamma)
  eigval[i,] <- ev$values
}

# plot eigenvalues of Gram matrices
all <- c(eigval[1,],eigval[2,],eigval[3,])
plot(1:n,all[1:n],col='royalblue',log='y',
     main='Sample Gram matrix eigenvalues',
     xlab='j',ylab='Eigenvalue',pch=20,
     ylim=range(all))
points(1:n,all[(1:n)+n],col='seagreen',pch=20)
points(1:n,all[(1:n)+2*n],col='salmon',pch=20)
print(paste("5-km condition number ",cond[1]))
print(paste("10-km condition number ",cond[2]))
print(paste("15-km condition number ",cond[3]))

# part (d)
hq <- (1:3)*5*1000
const <- (pi*(2*Gamma0*Dz)^2)/(2*hq)
diff <- outer(X, X, "-")
data1 <- gdata$Dataset1
N <- ((100+100)/dx)+1
xiq <- (((0:(N-1))*dx)-100)*1000
# Compute Gram matrix for each
Gamma5 <- const[1] / ((diff/(2*hq[1]))^2+1)
Gamma10 <- const[2] / ((diff/(2*hq[2]))^2+1)
Gamma15 <- const[3] / ((diff/(2*hq[3]))^2+1)
# Compute coefficients alpha for each, Model 1
alpha5 <- solve(Gamma5, data1)
alpha10 <- solve(Gamma10, data1)
alpha15 <- solve(Gamma15, data1)
# define data representers at each h as functions
g5 <- function(x,xi) {g(x,xi,hq[1],Dz)
}
g10 <- function(x,xi) {g(x,xi,hq[2],Dz)
}
g15 <- function(x,xi) {g(x,xi,hq[3],Dz)
}
# define matrices by which to multiply alpha to get drho
G5 <- t(outer(X,xiq,g5))
G10 <- t(outer(X,xiq,g10))
G15 <- t(outer(X,xiq,g15))
# compute drho using formula given
drho5 <- G5 %*% alpha5
drho10 <- G10 %*% alpha10
drho15 <- G15 %*% alpha15

# plot all
all <- c(drho5,drho10,drho15)
plot(xiq,drho5,type='l',col='goldenrod',xaxt="n",
     xlab='xi [km]', main='Best-fit density for Model 1',
     ylab='deltarho [kg^(1)m^(-3)]',ylim=range(all))
lines(xiq,drho10,'l',col='darkorange',xaxt="n")
lines(xiq,drho15,'l',col='deepskyblue',xaxt="n")
axis(1, at=((0:20)*20-100)*1000, labels=(0:20)*20-100)

# plot h=5 km case
plot(xiq,drho5,type='l',col='goldenrod',xaxt="n",
     xlab='xi [km]', main='Best-fit density for Model 1,
     h=5 km', ylab='deltarho [kg^(1)m^(-3)]')
axis(1, at=((0:20)*20-100)*1000, labels=(0:20)*20-100)

# plot h=10 km case
plot(xiq,drho10,type='l',col='darkorange',xaxt="n",
     xlab='xi [km]', main='Best-fit density for Model 1,
     h=10 km', ylab='deltarho [kg^(1)m^(-3)]')
axis(1, at=((0:20)*20-100)*1000, labels=(0:20)*20-100)

# plot h=15 km case
plot(xiq,drho15,type='l',col='deepskyblue',xaxt="n",
     xlab='xi [km]', main='Best-fit density for Model 1,
     h=15 km', ylab='deltarho [kg^(1)m^(-3)]')
axis(1, at=((0:20)*20-100)*1000, labels=(0:20)*20-100)

# part (e)
data2 <- gdata$Dataset2
# Compute coefficients alpha for each, Model 2
alpha5 <- solve(Gamma5, data2)
alpha10 <- solve(Gamma10, data2)
alpha15 <- solve(Gamma15, data2)
# recompute drho using formula given
drho5 <- G5 %*% alpha5
drho10 <- G10 %*% alpha10
drho15 <- G15 %*% alpha15

# plot all
all <- c(drho5,drho10,drho15)
plot(xiq,drho5,type='l',col='goldenrod',xaxt="n",
     xlab='xi [km]', main='Best-fit density for Model 2',
     ylab='deltarho [kg^(1)m^(-3)]',ylim=range(all))
lines(xiq,drho10,'l',col='darkorange',xaxt="n")
lines(xiq,drho15,'l',col='deepskyblue',xaxt="n")
axis(1, at=((0:20)*20-100)*1000, labels=(0:20)*20-100)

# plot h=5 km case
plot(xiq,drho5,type='l',col='goldenrod',xaxt="n",
     xlab='xi [km]', main='Best-fit density for Model 2,
     h=5 km', ylab='deltarho [kg^(1)m^(-3)]')
axis(1, at=((0:20)*20-100)*1000, labels=(0:20)*20-100)

# plot h=10 km case
plot(xiq,drho10,type='l',col='darkorange',xaxt="n",
     xlab='xi [km]', main='Best-fit density for Model 2,
     h=10 km', ylab='deltarho [kg^(1)m^(-3)]')
axis(1, at=((0:20)*20-100)*1000, labels=(0:20)*20-100)

# plot h=15 km case
plot(xiq,drho15,type='l',col='deepskyblue',xaxt="n",
     xlab='xi [km]', main='Best-fit density for Model 2,
     h=15 km', ylab='deltarho [kg^(1)m^(-3)]')
axis(1, at=((0:20)*20-100)*1000, labels=(0:20)*20-100)