import matplotlib.pylab as plt
import numpy as np
import os

with open(os.path.normpath(os.getcwd()) + '/data/checkpoint256.txt') as f:
    lines = f.readlines()
    t = [float(line.split(',')[0]) for line in lines]
    x = [float(line.split(',')[1]) for line in lines]
    y = [float(line.split(',')[2]) for line in lines]
    u = [float(line.split(',')[3]) for line in lines]

num_zones = int(np.sqrt(len(x)))

x_values = np.zeros(num_zones)
y_values = np.zeros(num_zones)
u_values = np.zeros((num_zones, num_zones))

for i in range(num_zones):
    x_values[i] = x[num_zones * i]
    y_values[i] = y[i]
    for j in range(num_zones):
        u_values[j, i]  = u[num_zones * i + j]

fig, ax = plt.subplots(1, 1, figsize=(6, 6))
ax.pcolormesh(x_values, y_values, u_values, shading='auto', cmap='plasma', vmin=0, vmax=4E4)
plt.tight_layout()
plt.show()
