# define setup code
setup = '''
# define parameters
N = 10000 # batch size
sigma = 1
mu = 0
a = -5 # start value for domain
b = 5 # end value for domain

# download normal distribution PDF from JAX
from jax.scipy.stats.norm import pdf

# define function f from Problem 5, evaluated at points in x
gaussian_f = lambda x: pdf(x, loc=mu, scale=sigma)[0]  # evaluate f at points in x

# define first-derivative f' of f from Problem 5, evaluated at the parameters
# define differentiation via jax.grad, do over multiple arguments using vmap
from jax import grad, vmap
gaussian_f_prime = vmap(grad(gaussian_f))

from jax.numpy import linspace
# use gaussian_f_prime to predict f'(x) for N uniform points in [a, b]
X = linspace(a, b, N).reshape((N, 1))
'''

# define test code
test = "gaussian_f_prime(X)"

# test_2 = "pdf(X, loc=mu, scale=sigma) * (mu - X) / (sigma ** 2)"

# test timing
n = int(5e3) # number of trials
import timeit # get timer to time n runs of gaussian_f 
total_time = timeit.timeit(setup=setup, stmt=test, number=n)
print(f'Average time: {total_time / n} seconds')

# import timeit # get timer to time n runs of gaussian_f 
# total_time = timeit.timeit(setup=setup, stmt=test_2, number=n)
# print(f'Average time: {total_time / n} seconds')

exec(setup) # execute setup code now within script
f_prime_jax = gaussian_f_prime(X)
# compute f'(x) using the exact formula
f_prime_exact = pdf(X, loc=mu, scale=sigma) * (mu - X) / (sigma ** 2)

# plot f_prime_jax and f_prime_exact
import matplotlib.pyplot as plt
plt.plot(X, f_prime_jax, label='Via AD')
plt.plot(X, f_prime_exact, label='Via exact formula')
plt.title("Comparison between results of AD and exact formula for $f'(x)$")
plt.xlabel('x')
plt.ylabel("$f'(x)$")
plt.legend()
plt.savefig('6_1.pdf')
plt.close()

from jax.numpy import abs
plt.plot(X, abs(f_prime_exact-f_prime_jax))
plt.title("Comparison between results of AD and exact formula for $f'(x)$")
plt.ylabel(r"$\left|f_{\mathrm{exact}}'(x)-f_{\mathrm{AD}}'(x)\right|$")
plt.xlabel('x')
plt.savefig('6_2.pdf')
plt.close()