#!/usr/bin/env -S python3 -i

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata


import os
cwd = os.path.dirname(__file__)
fn = f"{cwd}/my.mplstyle"
plt.style.use(fn)


print("Loading data...")
df = pd.read_csv('model/pressure.csv.gz', compression='gzip')

# Add colormap for 10th timestep as background
fig, ax = plt.subplots(figsize=(12, 8))

x_bg  = df['x'].values
z_bg  = df['z'].values
dt_bg = df['PRES'].values

x0 = x_bg[0]
df = df[df['x']==x0]

ax.scatter(df.PRES, df.z)

# # Create fine grid
# x_min, x_max = x_bg.min(), x_bg.max()
# z_min, z_max = z_bg.min(), z_bg.max()
# xi = np.logspace(np.log10(x_min), np.log10(x_max), 100)
# zi = np.linspace(z_bg.min(), z_bg.max(), 1000)
# X_bg, Z_bg = np.meshgrid(xi, zi)

# DT_bg = griddata((x_bg, z_bg), dt_bg, (X_bg, Z_bg), method='linear')
# contourf = ax.contourf(X_bg, Z_bg, DT_bg, levels=30, alpha=0.7, vmin=60e6)
# plt.colorbar(contourf, ax=ax, label=r'$\Delta P$')

# ax.set_xlim(0,350)
# ax.set_ylim(5450,5550)
ax.invert_yaxis()

plt.show()

