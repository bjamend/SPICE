import matplotlib.pylab as plt
import numpy as np
import os

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.sans-serif": ["Times New Roman"]})

def g(x):
    return 0.4 / x

with open(os.path.normpath(os.getcwd()) + '/data/128/checkpoint0.txt') as f:
    lines = f.readlines()
    u = [float(line.split(',')[3]) for line in lines]
    m = np.mean(u)
    for i in range(len(u)):
        u[i] = u[i] - m

with open(os.path.normpath(os.getcwd()) + '/data/256/checkpoint0.txt') as f:
    lines = f.readlines()
    u1 = [float(line.split(',')[3]) for line in lines]
    m1 = np.mean(u1)
    for i in range(len(u1)):
        u1[i] = u1[i] - m1

with open(os.path.normpath(os.getcwd()) + '/data/512/checkpoint0.txt') as f:
    lines = f.readlines()
    u2 = [float(line.split(',')[3]) for line in lines]
    m2 = np.mean(u2)
    for i in range(len(u2)):
        u2[i] = u2[i] - m2

with open(os.path.normpath(os.getcwd()) + '/data/1024/checkpoint0.txt') as f:
    lines = f.readlines()
    u3 = [float(line.split(',')[3]) for line in lines]
    m3 = np.mean(u3)
    for i in range(len(u3)):
        u3[i] = u3[i] - m3

print(m, m1, m2, m3)

bin_means = (np.histogram(u, bins=np.linspace(-30, 40.0, 500), density='True')[0])
bin_bounds = (np.histogram(u, bins=np.linspace(-30.0, 40.0, 500), density='True')[1])

bin_means1 = (np.histogram(u1, bins=np.linspace(-30, 40.0, 500), density='True')[0])
bin_bounds1 = (np.histogram(u1, bins=np.linspace(-30.0, 40.0, 500), density='True')[1])

bin_means2 = (np.histogram(u2, bins=np.linspace(-30, 40.0, 500), density='True')[0])
bin_bounds2 = (np.histogram(u2, bins=np.linspace(-30.0, 40.0, 500), density='True')[1])

bin_means3 = (np.histogram(u3, bins=np.linspace(-30, 40.0, 500), density='True')[0])
bin_bounds3 = (np.histogram(u3, bins=np.linspace(-30.0, 40.0, 500), density='True')[1])

plt.step(bin_bounds[:-1], bin_means, color='k', label='$n=128$')
plt.step(bin_bounds1[:-1], bin_means1, color='r', label='$n=256$')
plt.step(bin_bounds2[:-1], bin_means2, color='g', label='$n=512$')
plt.step(bin_bounds3[:-1], bin_means3, color='b', label='$n=1024$')
plt.yscale('log')
plt.legend()
plt.xlim(-30.0, 40.0)
plt.ylim(5E-5, 2E0)
plt.xlabel('$X-\overline{X}$')
plt.ylabel('$P(X-\overline{X})$')
plt.tight_layout()
plt.show()
