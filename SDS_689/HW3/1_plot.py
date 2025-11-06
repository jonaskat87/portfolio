import sys
import numpy as np

data = np.loadtxt(sys.stdin, dtype=np.float64)
epochs = np.arange(len(data)) + 1

import matplotlib.pyplot as plt
plt.plot(epochs,data[:,0])
plt.title('Cross-entropy loss vs. epoch')
plt.xlabel('Epoch')
plt.ylabel('Cross-entropy loss')
plt.yscale('log')
plt.savefig('1_loss.pdf')
plt.close()

plt.plot(epochs,data[:,1])
plt.title('Training error vs. epoch')
plt.xlabel('Epoch')
plt.ylabel('Training error (fraction of misidentified digits)')
plt.yscale('log')
plt.savefig('1_error.pdf')
plt.close()