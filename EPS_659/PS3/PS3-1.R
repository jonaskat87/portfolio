N <- 9
xq <- 1 : N
seq <- xq * xq

avg = sum(seq) / N # sample mean
med = median(seq) # sample median
hmean = 1 / (sum(1 / seq) / N) # sample harmonic mean
gmean = exp(sum(log(seq)) / N) # sample geometric mean

print(avg)
print(med)
print(hmean)
print(gmean)
