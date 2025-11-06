# the following function computes the L2 norm of a (1D) array
# does not use NumPy, but uses built-in math module in Python

import math

def L2_norm(v):
    norm_sq = [vi ** 2 for vi in v]
    return math.sqrt(sum(norm_sq))

# test L2_norm of provided vector
print(L2_norm([math.cos(math.pi/4)*math.sin(math.pi/8), 
               math.sin(math.pi/4)*math.sin(math.pi/8),
               math.cos(math.pi/8)]))