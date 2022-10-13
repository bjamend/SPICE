import matplotlib.pylab as plt
import numpy as np
import os, os.path

DIR = './data'
N = len([name for name in os.listdir(DIR) if os.path.isfile(os.path.join(DIR, name))]) - 1

for l in range(N):
    i = 10 * l
    f = open(f"./data/data{i}.txt",'r')
    x = []
    y = []
    u = []

    for row in f:
        row = row.split(',')
        x.append(float(row[0]))
        y.append(float(row[1]))
        u.append(float(row[2]))

        num_zones = int(np.sqrt(len(x)))

    x_values  = np.zeros(num_zones)
    y_values  = np.zeros(num_zones)
    u_values  = np.zeros((num_zones, num_zones))
    u0_values = np.zeros((num_zones, num_zones))

    for k in range(num_zones):
        x_values[k] = x[num_zones * k]
        y_values[k] = y[k]
        for j in range(num_zones):
            u_values[j, k]  = u[num_zones * k + j]

    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    ax.pcolormesh(x_values, y_values, u_values, shading='auto', cmap='plasma', vmin=0.0, vmax=0.15)
    fig.tight_layout()
    fig.savefig(f'./frames/{i}.png')

    f.close()
