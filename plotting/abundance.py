import numpy as np
import matplotlib
import matplotlib.pylab as plt
import random

matplotlib.use('Agg')

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.sans-serif": ["Times New Roman"]})

N = 164
num_zones = 128
H = 2.80E59
solar_fe_h  = -4.52
solar_eu_mg = -7.02

stars_fe_h  = []
stars_eu_mg = []

iteration = 0;

for i in range(N):
    f1 = open(f"./data/SNe data/data{1000*i}.txt", 'r')
    f2 = open(f"./data/BNS data/data{1000*i}.txt", 'r')
    Fe = []
    Mg = []
    Eu = []

    for row in f1:
        row = row.split(',')
        Fe.append(float(row[3]) * 4.76E55)
        Mg.append(float(row[3]) * 5.96E55)

    for row in f2:
        row = row.split(',')
        Eu.append(float(row[3]) * 3.57E49)

    count = 0;

    while (count < 4):

        for j in range(i):
            x_i = random.randint(0, num_zones - 1)
            y_i = random.randint(0, num_zones - 1)
            Fe_H = np.log10((Fe[num_zones * x_i + y_i] + 1.0E-10) / H)
            Eu_Mg = np.log10((Eu[num_zones * x_i + y_i] + 1.0E-10) / (Mg[num_zones * x_i + y_i] + 1.0E-10))
            stars_fe_h.append(Fe_H - solar_fe_h)
            stars_eu_mg.append(Eu_Mg - solar_eu_mg)

        plt.scatter(stars_fe_h, stars_eu_mg, color='k', s=5, alpha=1.0, edgecolors='none')
        plt.xlim(-2.5, 1.25)
        plt.ylim(-1.0, 1.75)
        plt.xlabel(r'[Fe/H]')
        plt.ylabel(r'[Eu/Mg]')
        plt.savefig(f"./frames/{iteration}.png")
        plt.close()

        count += 1
        iteration += 1

    print(i)
