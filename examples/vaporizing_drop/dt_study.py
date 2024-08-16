# Import libraries
import os
import numpy as np
import matplotlib.pyplot as plt
import math

# Reference values
D0    = 0.01
R0    = 0.5 * D0
mflux = 5
rhoL  = 1000
muL   = 1.05e-3
# tref  = rhoL * D0**2 / muL
tR  = R0 * rhoL / mflux

# Arrays
dt  = []
t   = []
r   = []
cfl = []

# Get the data
path = './monitor'
for filename in os.listdir(path):
    if 'simulation_dt_' in filename:
        data = np.loadtxt(path + '/' + filename, skiprows=2)
        t.append(data[:, 1 ] / tR)
        r.append(data[:, 14] / R0)
        cfl.append(data[:, 3].max())
        dt.append(float(filename.replace('simulation_dt_', '')))

# Sort
sorted_lists = sorted(zip(dt, t))
dummy, t = zip(*sorted_lists)
sorted_lists = sorted(zip(dt, r))
dummy, r = zip(*sorted_lists)
sorted_lists = sorted(zip(dt, cfl))
dt, cfl = zip(*sorted_lists)

# Analytical solution
t_ext = np.linspace(start=0, stop=1e-2, num=100)
r_ext = (R0 - mflux/rhoL * t_ext) / R0
t_ext = t_ext / tR

# Make colors
colors = plt.cm.plasma(np.linspace(0, 1, len(dt)))

# Use latex font
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

# Font sizes
fnt1 = 16
fnt2 = 14
fnt3 = 10

# Create a figure
fig, ax = plt.subplots(1, 1, figsize=(3.75, 3.5))

# Plot
plt.plot(t_ext, r_ext, ls='-', lw=2, color='k', label=r'$Analytical$')
for i, delta_t in enumerate(dt):
    plt.plot(t[i], r[i], ls='--', lw=2, color=colors[i], label=r'$\Delta t = \,$' + rf'${delta_t:.2e}$')# + r'$(s)$' + r'$, \, CFL_{max} = \,$' + rf'${cfl[i]:.2f}$')

plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel(r'$t/t_{ref}$', fontsize=fnt1)
plt.ylabel(r'$r/R_{0}$', fontsize=fnt1)
plt.xticks(fontsize=fnt2)
plt.yticks(fontsize=fnt2)
plt.grid(which='major', axis='both', color='k', linestyle='--', linewidth=0.6, alpha=0.25)
plt.legend(frameon=False, loc='upper right', fontsize=fnt3)
plt.tight_layout()
plt.savefig('./dt_study.pdf')