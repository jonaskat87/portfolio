# plots MSE loss vs. epoch for different widths

import sys
import numpy as np
data = np.genfromtxt(sys.stdin, dtype=np.float64, invalid_raise = False)
widths = data[-1] # network widths
errors = data[-2] # L2 relative errors for testing
losses = data[:-2]

import matplotlib.pyplot as plt
plt.title('Training across different network widths')
plt.xlabel('Epochs')
plt.ylabel('MSE loss (linear scale)')
for i, w in enumerate(widths):
    plt.plot(np.arange(len(losses)), losses[:,i],
             label='Width: '+str(w)+', Error: '+str(100*errors[i])+'%')
plt.legend()
plt.savefig('2_linear.pdf')
plt.close()

plt.title('Training across different network widths')
plt.xlabel('Epochs')
plt.ylabel('MSE loss (log scale)')
for i, w in enumerate(widths):
    plt.plot(np.arange(len(losses)), losses[:,i],
             label='Width: '+str(w)+', Error: '+str(100*errors[i])+'%')
plt.legend()
plt.yscale('log')
plt.savefig('2_log.pdf')
plt.close()