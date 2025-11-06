# script for plotting the relative errors 
# input is text file where each row is a network width
# each column corresponds to a new run, 
# except the first which contains the list of widths

import sys
import numpy as np

err_data = np.loadtxt(sys.stdin, dtype=np.float64) # load error data
# preallocate where to store widths, error means, and error standard deviations
err = np.zeros((len(err_data), 3)) 
for i, w in enumerate(err_data):
    errs = w[1:]
    err[i] = np.array([w[0], np.mean(errs), np.std(errs)])

# save two plots, first in linear scaling and second with log scaling for both axes
import matplotlib.pyplot as plt
plt.errorbar(x=err[:,0], y=err[:,1], yerr=err[:,2], fmt='b-', ecolor='cyan')
plt.title('$L^{2}$ relative error vs. network width')
plt.xlabel('Network width')
plt.ylabel('$L^{2}$ relative error (linear scale)')
plt.tight_layout()
plt.savefig('3_lin_scaling.pdf') # change file name here
plt.close()

plt.errorbar(x=err[:,0], y=err[:,1], yerr=err[:,2], fmt='b-', ecolor='cyan')
plt.title('$L^{2}$ relative error vs. network width')
plt.xlabel('Network width (log scale)')
plt.ylabel('$L^{2}$ relative error (log scale)')
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.savefig('3_log_scaling.pdf') # change file name here
plt.close()
