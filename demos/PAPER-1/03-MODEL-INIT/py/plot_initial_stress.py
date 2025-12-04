#!/usr/bin/env -S python3 -i

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import os
cwd = os.path.dirname(__file__)
fn = f"{cwd}/my.mplstyle"
plt.style.use(fn)

df = pd.read_pickle("pkl/stress_profile.pkl")

# Capture the discontinuity
z_target = 0
idx_below = df[df['Z'] <= z_target].index[-1]
idx_above = df[df['Z'] > z_target].index[0]
sigma_xx_below = df.loc[idx_below, 'sigma_xx']
sigma_xx_above = df.loc[idx_above, 'sigma_xx']
z_below = df.loc[idx_below, 'Z']
z_above = df.loc[idx_above, 'Z']
difference = sigma_xx_above - sigma_xx_below
mid_sigma = (sigma_xx_below + sigma_xx_above) / 2
print(f"Z below: {z_below} m, sigma_xx: {sigma_xx_below:.2e} Pa")
print(f"Z above: {z_above} m, sigma_xx: {sigma_xx_above:.2e} Pa")
print(f"Difference: {difference/1e6:.2e} MPa")

# Add a new value for pressure to avoid interpolation errors
p0 = df.loc[idx_above, 'Pressure']

new_row = {
    'Z': 0.0,
    'X': None,
    'Y': None,
    'sigma_xx': None,
    'sigma_yy': None,
    'sigma_zz': None,
    'Pressure': p0,  # Use pressure from Z > 0
    'Temperature': None,
    'Delta_P': None,
    'Delta_T': None,
}
df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)
df = df.sort_values('Z').reset_index(drop=True)



# Create the plot
fig, [ax,ax2] = plt.subplots(1,2,figsize=(4, 4),sharey=True)

# Plot each stress component
ax.plot(df['sigma_xx']/1e6, df['Z']-5500, label=r'$\sigma_{xx}$', linewidth=1, ls=':')
ax.plot(df['sigma_yy']/1e6, df['Z']-5500, label=r'$\sigma_{yy}$', linewidth=1, ls='--')
ax.plot(df['sigma_zz']/1e6, df['Z']-5500, label=r'$\sigma_{zz}$', linewidth=1)

ax.text(mid_sigma/1e6, z_target-5504, r'$\Delta \sigma$'+f' = {difference/1e6:.2f} MPa', 
        fontsize=10, ha='center', va='center',
        bbox=dict(facecolor='white', edgecolor='None', alpha=0.8))


# Set labels and limits
ax.set_xlabel(r'$\sigma$ (MPa)')
ax.set_ylabel('Depth (m)')
ax.set_ylim(-5450, -5550)
ax.grid(True)
ax.legend(fontsize=12, loc='upper left', bbox_to_anchor=(0.03,0.99))
ax.invert_xaxis()

ax2.plot(df['Pressure']/1e6, df['Z']-5500, label=r'Pressure', linewidth=1)
ax2.set_xlabel("Pressure (MPa)")

plt.show()
