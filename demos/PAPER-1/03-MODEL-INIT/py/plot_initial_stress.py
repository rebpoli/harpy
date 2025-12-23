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


# Create the plot
fig, [ax,ax3] = plt.subplots(1,2,figsize=(9/2.54, 8/2.54),sharey=True)


# Plot each stress component
ax.plot(-df['sigma_xx']/1e6, df['Z']+5500, label=r'$\sigma_{h}$', linewidth=0.7, ls='--', c='r')
# ax.plot(df['sigma_yy']/1e6, df['Z']+5500, label=r'$\sigma_{yy}$', linewidth=0.7, ls='--')
ax.plot(-df['sigma_zz']/1e6, df['Z']+5500, label=r'$\sigma_{\text{v}}$', linewidth=0.7, c='k')


ax.text(-mid_sigma/1e6-3, z_target+5499, r'$\Delta \sigma$'+f' = {difference/1e6:.1f} MPa', 
        fontsize=7, ha='center', va='bottom', c='r');
#         bbox=dict(facecolor='k', edgecolor='None', alpha=0.05))

# ax.annotate('',
#             xy=(-mid_sigma/1e6-5, z_target+5498),  # Arrow points TO this location
#             xytext=(-mid_sigma/1e6, z_target),           # Arrow starts FROM this location (text position)
#             arrowprops=dict(
#                 arrowstyle='->',
#                 connectionstyle='arc3,rad=0.3',  # rad controls curvature
#                 lw=2,
#                 color='black'
#             ))

ax.plot(df['Pressure']/1e6, df['Z']+5500, label=r'$p$', linewidth=0.7, ls='-.', c='b')


# Set labels and limits
ax.set_xlabel(r'$-$Stress,  Pressure (MPa)')
ax.set_ylabel('Depth (m)')
ax.set_ylim(+5450, +5550)

ax.grid(True)
ax.legend(fontsize=7, loc='upper center', bbox_to_anchor=(0.5,0.99))
# ax.invert_xaxis()
ax.invert_yaxis()


# ax.set_xticks([60.0, 60.1, 60.2, 60.3])
ax3.plot(df['Temperature']-273, df['Z']+5500, label=r'Temperature', linewidth=0.7, ls='--', c='green')

xlim=ax3.get_xlim()
xrng=xlim[1]-xlim[0]
x = xlim[1]-0.03*xrng
ax3.text(x, 5499,"Caprock",   va='bottom', ha='right', fontsize=7)
ax3.text(x, 5501,"Reservoir", va='top',    ha='right', fontsize=7)


ax3.set_xlabel(r"Temperature ($^\circ\!$C)")

# Plot the reservoir shades
for _ax in [ax, ax3] :
    _ax.axhspan(+5500, +5600, facecolor='k', alpha=0.2, zorder=0)


from matplotlib.ticker import FuncFormatter
def smart_format(x, pos):
    # Format with 2 decimals, then remove trailing zeros
    s = f'{x:.2f}'
    return s.rstrip('0').rstrip('.')
# for ax in [ax2]:
#     ax.xaxis.set_major_formatter(FuncFormatter(smart_format))

fig.savefig("png/init_stress_p_t.png", dpi=500)
#
# Manually set equal positions for all axes
# left_margin = 0.08
# right_margin = 0.95
# bottom = 0.15
# top = 0.90
# spacing = 0.02  # Space between subplots

# total_width = right_margin - left_margin - 2 * spacing
# subplot_width = total_width / 3
# ax.set_position([left_margin, bottom, subplot_width, top - bottom])
# ax2.set_position([left_margin + subplot_width + spacing, bottom, subplot_width, top - bottom])
# ax3.set_position([left_margin + 2 * (subplot_width + spacing), bottom, subplot_width, top - bottom])

# fig.subplots_adjust(wspace=0.1)  # Adjust horizontal spacing
# fig.tight_layout()
plt.show()
