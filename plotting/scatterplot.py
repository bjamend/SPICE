import matplotlib.pylab as plt
import random
import numpy as np

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.sans-serif": ["Times New Roman"]})

N = 262
num_zones = 256
H = 3E4
solar_fe_h  = -4.52

stars_fe_h  = []
times = []

for i in range(N):
    f = open(f"./data/checkpoint{i}.txt", 'r')
    temptimes = []
    ironvals  = []
    for row in f:
        row = row.split(',')
        temptimes.append(float(row[0]))
        ironvals.append(float(row[3]))
    for j in range(5):
        times.append(temptimes[0])
        x_index = random.randint(0, num_zones-1)
        y_index = random.randint(0, num_zones-1)
        stars_fe_h.append(np.log10(ironvals[x_index * num_zones + y_index] / 56.0 / H) - solar_fe_h)
    print(i)

plt.scatter(times[::-1], stars_fe_h, s=5, color='k')
plt.xlim(0.0, 10.0)
plt.ylim(-3.0, 1.0)
plt.xlabel('Age (Gyr)')
plt.ylabel('[Fe/H]')

plt.show()
