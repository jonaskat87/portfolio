# define setup code
setup = '''
# define parameters
N = 10000 # batch size
sigma = 1
mu = 0
a = -5 # start value for domain
b = 5 # end value for domain

# download normal distribution PDF and np.linspace from JAX
from jax.scipy.stats.norm import pdf
from jax.numpy import linspace

# define function f(x) from Problem 5
def gaussian_f(N, mu, sigma, a, b):
    # define grid on which to evaluate normal distribution, reshape to (N,1)
    X = linspace(a, b, N).reshape((N, 1)) 
    return pdf(X, loc=mu, scale=sigma)  # evaluate f(x) at points in X
'''

# define test code
test = "gaussian_f(N, mu, sigma, a, b)"

# test timing
n = int(5e4) # number of trials
import timeit # get timer to time n runs of gaussian_f 
total_time = timeit.timeit(setup=setup, stmt=test, number=n)
print(f'Average time: {total_time / n} seconds')
import platform # get CPU information
print(f'CPU: {platform.processor()}')