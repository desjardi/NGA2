# Import libraries
import numpy as np
import matplotlib.pyplot as plt
import math

# Reference values
mflux = np.array([0.05, 0.5, 5, 50, 500])
r0    = 0.005
rhoL  = 1000
muL   = 1.05e-3
D0    = 2.0 * r0
tref  = rhoL * D0**2 / muL

# Load file 0.05
data = np.loadtxt('./monitor/simulationmflux0.05', skiprows=2)

# Extract time and volume for 0.05
t0_05 = data[:, 1 ]
r0_05 = data[:, 14]

# Analytical solution for 0.05
mflux[0] = 0.05
r_ext0_05 = r0 - mflux[0]/rhoL * t0_05

# Load file 0.5
data = np.loadtxt('./monitor/simulationmflux0.5', skiprows=2)

# Extract time and volume for 0.5
t0_5 = data[:, 1 ]
r0_5 = data[:, 14]

# Analytical solution for 0.5
mflux[1] = 0.5
r_ext0_5 = r0 - mflux[1]/rhoL * t0_5

# Load file 5
data = np.loadtxt('./monitor/simulationmflux5', skiprows=2)

# Extract time and volume for 5
t5 = data[:, 1 ]
r5 = data[:, 14]

# Analytical solution for 5
mflux[1] = 0.5
r_ext5 = r0 - mflux[2]/rhoL * t5

# Load file 50
data = np.loadtxt('./monitor/simulationmflux50', skiprows=2)

# Extract time and volume for 50
t50 = data[:, 1 ]
r50 = data[:, 14]

# Analytical solution for 50
r_ext50 = r0 - mflux[3]/rhoL * t50

# Font sizes
fnt1 = 16
fnt2 = 14

# Use latex font
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

# Create a figure
fig, axs = plt.subplots(2, 2, figsize=(8, 8))

# Plot
axs[0, 0].plot(t0_05/tref, r_ext0_05/D0, ls='-'  , lw=2, color='k')
axs[0, 0].plot(t0_05/tref, r0_05    /D0, ls='--' , lw=2, color='b')
axs[0, 1].plot(t0_5 /tref, r_ext0_5 /D0, ls='-'  , lw=2, color='k')
axs[0, 1].plot(t0_5 /tref, r0_5     /D0, ls='--' , lw=2, color='b')
axs[1, 0].plot(t5   /tref, r_ext5   /D0, ls='-'  , lw=2, color='k')
axs[1, 0].plot(t5   /tref, r5       /D0, ls='--' , lw=2, color='b')
axs[1, 1].plot(t50  /tref, r_ext50  /D0, ls='-'  , lw=2, color='k')
axs[1, 1].plot(t50  /tref, r50      /D0, ls='--' , lw=2, color='b')
for i, ax in enumerate(axs.flat):
    ax.set_ylim(0, 0.52)
    ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax.set_xlabel(r'$t/t_{ref}$', fontsize=fnt1)
    ax.set_ylabel(r'$r/D_{0}$', fontsize=fnt1)
    ax.tick_params(axis='x', labelsize=fnt2)
    ax.tick_params(axis='y', labelsize=fnt2)
    ax.grid(which='major', axis='both', color='k', linestyle='--', linewidth=0.6, alpha=0.25)
    value = mflux[i]
    text = r'$\dot{m}^{\prime \prime}=\,$' + rf'${value:.2f}$'
    ax.text(0.95, 0.95, text, 
           transform=ax.transAxes,
           fontsize=fnt2, 
           verticalalignment='top', 
           horizontalalignment='right')
for ax in [axs[0, 1], axs[1, 1]]:
    ax.set_ylabel('')
    ax.set_yticklabels([])
for ax in [axs[0, 0], axs[0, 1]]:
    ax.set_xlabel('')
legends = axs[0, 0].legend([r'$Analytical$', r'$Numerical$'], frameon=False, loc='lower left', fontsize=fnt2)
axs[0, 0].add_artist(legends)
plt.tight_layout()
plt.savefig('./mflux_study.pdf')