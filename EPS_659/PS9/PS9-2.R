# part (a)
time <- 0:99
xdata <- sin(2*pi*time/20)
spec <- fft(xdata,inverse=TRUE)/length(xdata)
pspec <- abs(spec)
respec <- Re(spec) 
imspec <- Im(spec)
freq <- (0:99)/100
plot(freq,pspec,col='purple',pch=20,xlab='f_m',
     ylab='|Y_m|',main='Absolute value of spectrum, part (a)')
plot(freq,respec,col='blue',pch=20,xlab='f_m',
     ylab='Re{Y_m}',main='Real part of spectrum, part (a)')
plot(freq,imspec,col='red',pch=20,xlab='f_m',
     ylab='Im{Y_m}',main='Imaginary part of spectrum, part (a)')

# part (b)
xdata <- sin(2*pi*time/7)
spec <- fft(xdata,inverse=TRUE)/length(xdata)
pspec <- abs(spec)
respec <- Re(spec) 
imspec <- Im(spec)
freq <- (0:99)/100
plot(freq,pspec,col='purple',pch=20,xlab='f_m',
     ylab='|Y_m|',main='Absolute value of spectrum, part (b)')
plot(freq,respec,col='blue',pch=20,xlab='f_m',
     ylab='Re{Y_m}',main='Real part of spectrum, part (b)')
plot(freq,imspec,col='red',pch=20,xlab='f_m',
     ylab='Im{Y_m}',main='Imaginary part of spectrum, part (b)')

# part (c)
xdata <- sin(2*pi*time/20)
xpad <- rep(0.0,1000)
xpad[1:100] <- xdata
specpad <- fft(xpad,inverse=TRUE)/length(xdata)
fpad <- (0:999)/1000
plot(fpad,abs(specpad),col='darkorchid',pch=20,
     main='Absolute value of spectrum, part (c)')
plot(fpad,Re(specpad),col='deepskyblue',pch=20,
     main='Real part of spectrum, part (c)')
plot(fpad,Im(specpad),col='orangered',pch=20,
     main='Imaginary part of spectrum, part (c)')

plot(fpad[1:100],abs(specpad)[1:100],col='darkorchid',
     type="l",main='Spectrum (all components) for part (c),
      0<=fpad<0.1',xlab='fpad (between 0 and 0.1)',
     ylim=range(c(abs(specpad)[1:100],Re(specpad)[1:100],
                  Im(specpad)[1:100])),
     ylab='Spectrum (abs(specpad), Re(specpad), or Im(specpad)')
lines(fpad[1:100],Re(specpad)[1:100],col='deepskyblue')
lines(fpad[1:100],Im(specpad)[1:100],col='orangered')