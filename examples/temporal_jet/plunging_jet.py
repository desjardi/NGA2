import pandas as pd
from io import StringIO
import matplotlib.pyplot as plt
import numpy as np

# --- Raw tabulated data ---
raw_data = """\
1      2.360e-2   5.0      2.00    0
2      2.360e-2   7.5      9.00    0
3      8.290e-2   12.5     2.00    0.5
4      8.290e-2   2.5      3.25    0
5      8.290e-2   2.5      5.00    0
6      8.290e-2   20.0     5.00    1
7      8.290e-2   2.5      9.00    0
8      8.290e-2   7.5      9.00    1
9      1.623e-1   2.5      9.00    0.5
10     1.623e-1   5.0      9.00    1
11     1.623e-1   7.5      9.00    1
12     2.130e-1   1.8      3.25    0
13     2.130e-1   2.5      9.00    1
14     2.130e-1   5.0      8.00    1
15     1.623e-1   2.5      8.00    0
"""

# --- Load into a DataFrame (no header) ---
df = pd.read_csv(StringIO(raw_data), sep=r'\s+', header=None)
df.columns = ['Case', 'D', 'U', 'H', 'Aerated']
df = df.astype({'D': float, 'U': float, 'H': float, 'Aerated': float})

# --- Fluid properties (water) ---
rho = 1000.0       # kg/m^3
mu = 1.137e-3      # Pa*s (dynamic viscosity)
sigma = 0.0728     # N/m (surface tension)
g = 9.81           # m/s^2

# --- Dimensionless numbers ---
df['Re'] = rho * df['U'] * df['D'] / mu
df['We'] = rho * df['U']**2 * df['D'] / sigma
df['Fr'] = df['U'] / np.sqrt(g * df['D'])
df['Oh'] = mu / np.sqrt(rho * sigma * df['D'])

# --- Marker styles ---
marker_map = {
    0.0: ('o', 'Non-aerated'),
    0.5: ('s', 'Weakly aerated'),
    1.0: ('^', 'Aerated')
}

# --- Scatter Plot (log-log) ---
plt.figure(figsize=(8, 6))

for value, (marker, label) in marker_map.items():
    subset = df[df['Aerated'] == value]
    plt.scatter(
        subset['Re'], subset['We'],
        label=label,
        marker=marker,
        s=80,
        edgecolor='black',
        alpha=0.85
    )

# --- Fr = Fr_value line ---
Fr_value = 10  # Change this value as needed
D_vals = np.logspace(np.log10(df['D'].min()*0.5),
                     np.log10(df['D'].max()*2), 300)
Re_line = rho * D_vals * Fr_value * np.sqrt(g * D_vals) / mu
We_line = rho * (Fr_value**2) * g * D_vals**2 / sigma
plt.plot(Re_line, We_line, 'k--', label=f'Fr = {Fr_value}')

# --- Log scales ---
plt.xscale('log')
plt.yscale('log')

# --- Tight axis limits ---
x_min = df['Re'].min() * 0.7
x_max = df['Re'].max() * 1.3
y_min = df['We'].min() * 0.7
y_max = df['We'].max() * 1.3
plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)

# --- Labels and legend ---
plt.xlabel('Reynolds number (Re)')
plt.ylabel('Weber number (We)')
plt.title(f'Weber vs Reynolds (log-log) with Fr = {Fr_value} line')
plt.legend()
plt.grid(True, which='both', linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()