import matplotlib
import matplotlib.pylab as plt
import numpy as np
import os, os.path

# matplotlib.use('Agg')

DIR = './data'
N = len([name for name in os.listdir(DIR) if os.path.isfile(os.path.join(DIR, name))])

means = []
counters = []

counter = 0

for i in range(N):
    with open(os.path.normpath(os.getcwd()) + f'/data/checkpoint{i}.txt') as f:
        lines = f.readlines()
        u = [float(line.split(',')[3]) for line in lines]

    m = np.mean(u)
    print(m)
    means.append(m)
    counters.append(counter)
    counter += 1

    # fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    # plt.hist(u, bins=np.logspace(np.log10(0.001),np.log10(10.0), 50), density='True', color='k')
    # ax.set_yscale('log')
    # ax.set_xscale('log')
    # ax.set_xlim(1E-3, 1E1)
    # ax.set_ylim(1E-2, 1E2)
    # plt.tight_layout()
    # fig.savefig(f'./frames/{i}.png')
    f.close()

plt.scatter(counters, means)
plt.show()
