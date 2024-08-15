# Import libraries
import numpy as np
import matplotlib.pyplot as plt
import math

# Reference values
D0    = 0.01
R0    = D0 / 2
mflux = 5
rhoL  = 1000
muL   = 1.05e-3
# tref  = rhoL * D0**2 / muL
tR  = R0 * rhoL / mflux

# Load file
data = np.loadtxt('./monitor/simulation', skiprows=2)

# Extract time and volume
t = data[:, 1 ]
r = data[:, 14]

# Analytical solution
r_ext = R0 - mflux/rhoL * t

# Use latex font
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

# Create a figure
fig, ax = plt.subplots(1, 1, figsize=(4, 4))

# Plot
plt.plot(t/tR, r_ext/R0, ls='-'  , lw=2, color='k')
plt.plot(t/tR, r    /R0, ls='--' , lw=2, color='b')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel(r'$t/t_{ref}$', fontsize=12)
plt.ylabel(r'$r/R_{0}$', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(which='major', axis='both', color='k', linestyle='--', linewidth=0.6, alpha=0.25)
legends = ax.legend([r'$Analytical$', r'$Numerical$'], frameon=False, loc='lower left', fontsize=10)
ax.add_artist(legends)
plt.tight_layout()
plt.savefig('./r_vs_t.pdf')