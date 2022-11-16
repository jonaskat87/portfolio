# Start the clock!
ptm <- proc.time()

n <- 4000 # number of realizations (set each time)

# preallocate error array
errors <- rep(0, 8) 
# preallocate statistic values (exact and empirical)
stat <- matrix(0, nrow=2, ncol=8)
m <- seq(-5, 5, by=1) # sample values of mean
s <- seq(1, 5, by=1) # sample values of standard deviation
for (mx in m) {
  for (my in m) {
    for (sx in s) {
      for (sy in s) {
        stat[1, 1] <- mx
        stat[1, 2] <- mx^2 + sx^2
        stat[1, 3] <- my
        stat[1, 4] <- my^2 + sy^2
        stat[1, 5] <- 0.5 * (mx + my)
        stat[1, 6] <- 0.25 * (mx^2 + my^2 + sx^2 + sy^2) + 
          0.5*mx*my
        stat[1, 7] <- mx^2 + my^2 + sx^2 + sy^2 + 2*mx*my
        stat[1, 8] <- 2.25 * (mx^2 + sx^2) + 1.5*mx*my + 
          0.25*(my^2 + sy^2)
        x <- rnorm(n, mean=mx, sd=sx)
        y <- rnorm(n, mean=my, sd=sy)
        z <- (x+y)/2
        stat[2, 1] <- mean(x)
        stat[2, 2] <- mean(x^2)
        stat[2, 3] <- mean(y)
        stat[2, 4] <- mean(y^2)
        stat[2, 5] <- mean(z)
        stat[2, 6] <- mean(z^2)
        stat[2, 7] <- mean((x+y)^2)
        stat[2, 8] <- mean((x+z)^2)
        for (i in seq(1, 8, by=1)) {
          errors[i] <- errors[i] + (stat[1, i] - stat[2, i])^2
        }
      }
    }
  }
}
errors <- errors / ((length(m)*length(s))^2)
print(errors)

# Stop the clock
print(proc.time() - ptm)

