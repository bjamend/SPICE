import matplotlib.pylab as plt
import numpy as np
import os

def g(x):
    return 0.4 / x

x_vals = np.logspace(-2, 1, 50)

with open(os.path.normpath(os.getcwd()) + '/data/checkpoint122.txt') as f:
    lines = f.readlines()
    u = [float(line.split(',')[3]) for line in lines]

fig, ax = plt.subplots(1, 1, figsize=(6, 6))
plt.hist(u, bins=np.logspace(np.log10(0.01),np.log10(1.0), 50), density='True', color='k')
plt.plot(x_vals, g(x_vals))
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlim(1E-2, 1E0)
ax.set_ylim(1E-1, 1E2)
plt.tight_layout()
plt.show()
