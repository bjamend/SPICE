import matplotlib.pylab as plt
import numpy as np
import os

def g(x):
    return 0.4 / x

with open(os.path.normpath(os.getcwd()) + '/data/checkpoint200.txt') as f:
    lines = f.readlines()
    u = [float(line.split(',')[3]) for line in lines]
    m = np.mean(u)
    for i in range(len(u)):
        u[i] = u[i] - m

# bin_means = (np.histogram(u, bins=np.logspace(np.log10(0.0001),np.log10(1000.0), 250), density='True')[0])
# bin_bounds = (np.histogram(u, bins=np.logspace(np.log10(0.0001),np.log10(1000.0), 250), density='True')[1])

bin_means = (np.histogram(u, bins=np.linspace(-10, 50.0, 500))[0])
bin_bounds = (np.histogram(u, bins=np.linspace(-10.0, 50.0, 500))[1])

plt.step(bin_bounds[:-1], bin_means, color='k')
# plt.xscale('log')
plt.yscale('log')
plt.xlim(-8.0, 40.0)
plt.ylim(6E-1, 1E5)
plt.show()
