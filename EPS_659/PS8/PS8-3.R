# part (a)
file <- paste0("C:/Users/jonas/Documents",
               "/EPS 659/PS8/gravity_pset8.csv")
gdata <- read.csv(file)
X <- gdata$X
n <- length(X)

dx <- 0.01 # grid size for plot
h <- 5*1000 # anomaly depth
Dz <- 1000 # delta z
xi0q <- c(0,20,50)*1000 # different cases for xi
Gamma0 <- 6.67408e-11
# define data representer via a function
g <- function(x,xi,y,z) {2*Gamma0*y*z/(((x-xi)^2+y^2))
}
beta0 <- g(X,xi0q[1],h,Dz)
beta20 <- g(X,xi0q[2],h,Dz)
beta50 <- g(X,xi0q[3],h,Dz)
# compute Gram matrix Gamma
const <- (pi*(2*Gamma0*Dz)^2)/(2*h)
diff <- outer(X, X, "-") 
Gamma <- const / ((diff/(2*h))^2+1)
# solve for the coefficient vectors gamma
gamma0 <- solve(Gamma,beta0)
gamma20 <- solve(Gamma,beta20)
gamma50 <- solve(Gamma,beta50)

# sample values of xi for plotting
N <- ((100+100)/dx)+1
xiq <- (((0:(N-1))*dx)-100)*1000
# define data representer as a function
gmatrix <- function(x,xi) {g(x,xi,h,Dz)
}
# define matrices by which to multiply gamma
# to get resolution kernel
G <- t(outer(X,xiq,gmatrix))
# compute resolution kernel in each case for xi0
tdel0 <- G %*% gamma0
tdel20 <- G %*% gamma20
tdel50 <- G %*% gamma50

# plot all
all <- 1000*c(tdel0,tdel20,tdel50)
plot(xiq,1000*tdel0,type='l',col='indianred',xaxt="n",
     xlab='xi [km]',ylim=range(all),
     main='Sample resolution kernels for h=5 km',
     ylab='tildelta [km^(-1)]')
lines(xiq,1000*tdel20,'l',col='steelblue',xaxt="n")
lines(xiq,1000*tdel50,'l',col='goldenrod',xaxt="n")
axis(1, at=((0:10)*20-100)*1000, labels=(0:10)*20-100)

# part (b)
h <- 10*1000 # redefine anomaly depth at 10 km
# redefine all variables from part (a)
beta0 <- g(X,xi0q[1],h,Dz)
beta20 <- g(X,xi0q[2],h,Dz)
beta50 <- g(X,xi0q[3],h,Dz)
const <- (pi*(2*Gamma0*Dz)^2)/(2*h)
Gamma <- const / ((diff/(2*h))^2+1)
gamma0 <- solve(Gamma,beta0)
gamma20 <- solve(Gamma,beta20)
gamma50 <- solve(Gamma,beta50)
gmatrix <- function(x,xi) {g(x,xi,h,Dz)
}
G <- t(outer(X,xiq,gmatrix))
tdel0 <- G %*% gamma0
tdel20 <- G %*% gamma20
tdel50 <- G %*% gamma50

# plot all
all <- 1000*c(tdel0,tdel20,tdel50)
plot(xiq,1000*tdel0,type='l',col='indianred',xaxt="n",
     xlab='xi [km]',ylim=range(all),
     main='Sample resolution kernels for h=10 km',
     ylab='tildelta [km^(-1)]')
lines(xiq,1000*tdel20,'l',col='steelblue',xaxt="n")
lines(xiq,1000*tdel50,'l',col='goldenrod',xaxt="n")
axis(1, at=((0:10)*20-100)*1000, labels=(0:10)*20-100)