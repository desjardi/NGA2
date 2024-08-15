# Import libraries
import numpy as np
import matplotlib.pyplot as plt
import math

# Reference values
r0    = 0.005
mflux = 5
rhoL  = 1000
muL   = 1.05e-3
D0    = 2.0 * r0
tref  = rhoL * D0**2 / muL

# Load file 50
data = np.loadtxt('./monitor/simulation50', skiprows=2)

# Extract time and volume 50
t50 = data[:, 1 ]
r50 = data[:, 14]

# Load file 100
data = np.loadtxt('./monitor/simulation100', skiprows=2)

# Extract time and volume 100
t100 = data[:, 1 ]
r100 = data[:, 14]

# Load file 200
data = np.loadtxt('./monitor/simulation200', skiprows=2)

# Extract time and volume 200
t200 = data[:, 1 ]
r200 = data[:, 14]

# Load file 400
data = np.loadtxt('./monitor/simulation400', skiprows=2)

# Extract time and volume 400
t400 = data[:, 1 ]
r400 = data[:, 14]

# Analytical solution
r_ext = r0 - mflux * t50
# V_ext = math.pi * r_ext**2

# Make colors
colors = plt.cm.plasma(np.linspace(0, 1, 4))

# Use latex font
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

# Font sizes
fnt1 = 16
fnt2 = 14
fnt3 = 10

# Create a figure
fig, ax = plt.subplots(1, 1, figsize=(4, 4))

# Plot
plt.plot(t50/tref , r_ext /D0, ls='-'  , lw=2, color='k')
plt.plot(t50/tref , r50   /D0, ls='--' , lw=2, color=colors[0])
plt.plot(t100/tref, r100  /D0, ls='--' , lw=2, color=colors[1])
plt.plot(t200/tref, r200  /D0, ls='--' , lw=2, color=colors[2])
plt.plot(t400/tref, r400  /D0, ls='--' , lw=2, color=colors[3])
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel(r'$t/t_{ref}$', fontsize=fnt1)
plt.ylabel(r'$r/D_{0}$', fontsize=fnt1)
plt.xticks(fontsize=fnt2)
plt.yticks(fontsize=fnt2)
plt.grid(which='major', axis='both', color='k', linestyle='--', linewidth=0.6, alpha=0.25)
legends = ax.legend([r'$Analytical$', r'$50 \times 50$', r'$100 \times 100$', r'$200 \times 200$', r'$400 \times 400$'], frameon=False, loc='upper right', fontsize=fnt3)
ax.add_artist(legends)
plt.tight_layout()
plt.savefig('./mesh_study.pdf')