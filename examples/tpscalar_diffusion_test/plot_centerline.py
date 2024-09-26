# import libraries
import re
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import sys
sys.path.append('..')

# Function to read dat file
def file_to_data(file, data_names):

    # Read lines
    with open(file, 'r') as sf:
        lines = sf.readlines()

    # Store header
    names = lines[0].split()

    # Remove header line
    lines.pop(0)

    fdata = []
    for line in lines:
        split_line = line.split()
        fdata.append([float(val) for val in split_line])
    fdata = np.array(fdata)

    data_num = {}
    for name in data_names:
        data_num[name] = fdata[:, names.index(name)]
    
    return data_num

# Get the data
path      = './centerline'
arr_names = np.array(['x', 'Zl'])
data_num  = []
data_ext  = []
time_num  = []
time_ext  = []
n_times   = 0

pattern = r'[-+]?\d*\.\d+|\d+'
for filename in os.listdir(path):
    if filename.endswith('.dat'):
        match = re.search(pattern, filename)
        if 'numerical' in filename:
            n_times = n_times + 1
            time_num.append(float(match.group()))
            data_num.append(file_to_data(path +'/' + filename, arr_names))
        elif 'exact' in filename:
            time_ext.append(float(match.group()))
            data_ext.append(file_to_data(path + '/' + filename, arr_names))
        else:
            pass

sorted_lists = sorted(zip(time_num, data_num))
time_num, data_num = zip(*sorted_lists)

sorted_lists = sorted(zip(time_ext, data_ext))
time_ext, data_ext = zip(*sorted_lists)

# Use latex font
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

# Make colors
colors = plt.cm.plasma(np.linspace(0, 1, n_times))

# Create a figure
fig, ax = plt.subplots(1, 1, figsize=(3.0, 2.5))

# Loop over the times
for t_ind in range(n_times):
    ax.plot(data_ext[t_ind]['x'], 
            data_ext[t_ind]['Zl'], 
            color=colors[t_ind], 
            ls='-',
            lw=1.5,
            label = r'$t = {:.2f}\, [s]$'.format(time_ext[t_ind]))
    ax.plot(data_num[t_ind]['x'], 
            data_num[t_ind]['Zl'], 
            ls='none', 
            markerfacecolor=colors[t_ind],
            markeredgecolor=colors[t_ind],
            markeredgewidth=1,
            marker='o',
            ms=4,
            markevery=1)

plt.xlabel(r'$r$', fontsize=12)
plt.ylabel(r'$z_{l}$', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(which='major', axis='both', color='k', linestyle='--', linewidth=0.5, alpha=0.15)
# sm = plt.cm.ScalarMappable(cmap=plt.cm.plasma, norm=plt.Normalize(vmin=min(time_ext), vmax=max(time_ext)))
# sm.set_array([])
# cbar = plt.colorbar(sm, ax=ax, label=r'$t\, (s)$')
# cbar.set_ticks(time_ext)
# cbar.ax.tick_params(labelsize=10)
# cbar.set_label(r'$t\, (s)$', fontsize=12)
mk_lg = [Line2D([0], [0], color='k', linestyle='-', linewidth=1.5), 
         Line2D([0], [0], markerfacecolor='k', markeredgecolor='k', markersize=4, linestyle='None', marker='o', ms=4)]
legends = ax.legend(mk_lg, [r'$Exact$', r'$Numerical$'], frameon=False, loc='upper right', fontsize=10)
ax.add_artist(legends)
ax.legend(frameon=False, loc='center right', fontsize=10)
for spine in plt.gca().spines.values():
    spine.set_linewidth(0.75)
plt.tight_layout()
plt.xlim(0, 0.25)
plt.savefig('./centerline.pdf')