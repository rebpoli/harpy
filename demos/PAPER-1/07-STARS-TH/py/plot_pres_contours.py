#!/usr/bin/env -S python3 -i

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

from util import format_label
# from noflow_ghosts import add_noflow_ghosts, get_noflow_ghosts

import os
cwd = os.path.dirname(__file__)
fn = f"{cwd}/my.mplstyle"
plt.style.use(fn)

# ============================================================================

# Option 4: Auto-select evenly spaced timesteps (default)
AUTO_SELECT = 15  
contour_levels = [5e6]

# ============================================================================

# Load the gzipped CSV file
print("Loading data...")
df = pd.read_csv('pressure.csv.gz', compression='gzip')

# Get initial pressure field at t=0
df_t0 = df[df['t'] == 0].copy()
df_t0 = df_t0.rename(columns={'PRES': 'PRES_0'})
df_t0 = df_t0[['x', 'z', 'PRES_0']]

# We do not want very short noisy times
df = df[df.t > 0.1]

# Merge to get PRES_0 for all timesteps
df = df.merge(df_t0, on=['x', 'z'], how='left')

# Calculate DELTA_PRES
df['DELTA_PRES'] = df['PRES'] - df['PRES_0']

# Get unique timesteps (excluding t=0 since DELTA_PRES would be 0)
timesteps = sorted(df['t'].unique())
timesteps = [t for t in timesteps if t > 0]

print(f"\nTotal available timesteps: {len(timesteps)}")
print(f"Time range: {timesteps[0]:.2f} to {timesteps[-1]:.2e} days")


# Create logarithmically spaced time values
# t_min = 1*60*60#timesteps[0]
# t_max = 5*365*24*60*60#timesteps[-1]
t_min = 10*24*60*60
t_max = 10*365*24*60*60
log_times = np.logspace(np.log10(t_min), np.log10(t_max), AUTO_SELECT)
# log_times = np.linspace(t_min, t_max, AUTO_SELECT)

# log_times = []
# for t_h in [6, 12] :
#     log_times.append(t_h*60*60)
# for t_d in [1, 5, 10] :
#     log_times.append(t_d*24*60*60)
# for t_m in [1, 6, 11, 5*12, 10*12] :
#     log_times.append(t_m*30*24*60*60)

selected_timesteps = []
for log_t in log_times:
    idx = np.argmin(np.abs(np.array(timesteps) - log_t))
    selected_timesteps.append(timesteps[idx])
selected_timesteps = sorted(set(selected_timesteps))

sel_strs = [ format_label(t) for t in selected_timesteps ] 

print(f"Plotting timesteps: {sel_strs}")
print(f"Contour levels: {contour_levels}")

# Create the plot
fig, ax = plt.subplots(figsize=(9/2.54, 6/2.54))

# Define colormap for different timesteps
colors = plt.cm.viridis(np.linspace(0, 1, len(selected_timesteps)))

lab_list = [
    [ '10 days', (30, 5530), -88, 7 ],   #10 d
    [ ], # 15 d
    [ ], # 23 d
    [ ], # 1 m
    [ ], # 2 m
    [ ], # 2.6 m
    [ '4 months', (100, 5540), -80 , -5], # 4 m
    [ ], # 6.5 m
    [ ], # 9.5 m
    [ '15 months', (100, 5540), -72 , -5], # 15 m
    [ '23 months', (100, 5540), -70, -5 ], # 23 m
    [ '3 years', (300, 5530), -67 , -5], # 3 y
    [ ], # 4 y
    [ ], # 6.6 y
    [ ], # 10 y
]


# Plot contours for each selected timestep
for idx, t in enumerate(selected_timesteps):
    df_t = df[df['t'] == t]

#     x, z, delta_pres = add_noflow_ghosts( df_t, [5495, 5595] )
    
    # Get data
    x = df_t['x'].values
    z = df_t['z'].values
    delta_pres = df_t['DELTA_PRES'].values
    
    print(f"\nTimestep {t:.2e}: DELTA_PRES range = [{delta_pres.min():.2f}, {delta_pres.max():.2f}]")
    
    # Create a fine interpolation grid for smooth contours
    x_min, x_max = x.min(), x.max()
    z_min, z_max = z.min(), z.max()
    
    # Create a fine grid (increase resolution for smoother contours)
    n_points = 1000  # Number of points in each direction
    xi = np.logspace(np.log10(x_min), np.log10(x_max), n_points)
    zi = np.linspace(z_min, z_max, 2000)
    X, Z = np.meshgrid(xi, zi)
    
    # Interpolate DELTA_PRES onto the fine grid
    DELTA_PRES_grid = griddata((x, z), delta_pres, (X, Z), method='linear')
    
    # Check which contour levels are in range for this timestep
    valid_levels = [level for level in contour_levels 
                    if delta_pres.min() <= level <= delta_pres.max()]
    
    contour = ax.contour(X, Z, DELTA_PRES_grid, levels=valid_levels, 
                         colors='k', linestyles='--', linewidths=0.6)

    # Add inline labels to the contour
    l = lab_list[idx]
    if len(l) :
        labs = ax.clabel(contour, inline=False, fontsize=6, 
                          manual=[l[1]], fmt=l[0], inline_spacing=0,
                          )#, 
        labs[0].set_rotation(l[2])
        x,y = labs[0].get_position()
        theta = np.deg2rad(l[2])
        nx=np.sin(theta)
        ny=np.cos(theta)
        offset=l[3]

        labs[0].set_position((x+offset*nx,y+offset*ny))

ax.axhspan(5500, 5600, facecolor='gray', alpha=0.3, zorder=0)

# Labels and formatting
ax.set_xlabel('Distance from well (m)')
ax.set_ylabel('Depth (m)')
ax.set_xlim(0,350)
ax.set_ylim(5400,5650)
ax.invert_yaxis()

# Create title based on contour levels
# if len(contour_levels) == 1:
#     title = f'DELTA_PRES = {contour_levels[0]} Contours'
# else:
#     title = f'DELTA_PRES Contours at levels: {contour_levels}'

ax.set_title(r"Isobaric contours $\Delta P(t)=5$MPa, at selected times")

# Only show legend if we have contours
if ax.get_legend_handles_labels()[0]:
    ax.legend(loc='best')
ax.grid(True, alpha=0.3)


fig.savefig('png/pres_contours.png', dpi=300, bbox_inches='tight')


# # Add colormap for 10th timestep as background
# fig, ax = plt.subplots(figsize=(12, 8))
# if len(selected_timesteps) >= 10:
#     t_bg = selected_timesteps[5]  # 10th timestep (0-indexed)
#     print(f"T_BG: {t_bg}")
#     df_bg = df[df['t'] == t_bg]
#     
# #     x_bg, z_bg, dt_bg = add_noflow_ghosts( df_bg, [5495,5595])
#     x_bg  = df_bg['x'].values
#     z_bg  = df_bg['z'].values
#     dt_bg = df_bg['DELTA_PRES'].values

#     # Create fine grid
# #     xi = np.linspace(x_bg.min(), x_bg.max(), 1000)
#     x_min, x_max = x_bg.min(), x_bg.max()
#     z_min, z_max = z_bg.min(), z_bg.max()
#     xi = np.logspace(np.log10(x_min), np.log10(x_max), 100)
#     zi = np.linspace(z_bg.min(), z_bg.max(), 1000)
#     X_bg, Z_bg = np.meshgrid(xi, zi)

#     DT_bg = griddata((x_bg, z_bg), dt_bg, (X_bg, Z_bg), method='linear')
#     contourf = ax.contourf(X_bg, Z_bg, DT_bg, levels=30, cmap='RdBu_r', alpha=0.7)
#     plt.colorbar(contourf, ax=ax, label=r'$\Delta P$')

# #     contour = ax.contour(X_bg, Z_bg, DT_bg, levels=[2e6], colors='k', ls='--', linewidths=0.6)

# #     xg,zg,dg = get_noflow_ghosts( df_bg, [5495] )
# #     ax.scatter( xg, zg, c=dg, cmap='RdBu_r', vmin=0, vmax=5e6 )
# #     ax.scatter( x_bg, z_bg, marker='s', s=30, c=dt_bg, cmap='RdBu_r', vmin=0, vmax=5e6 )
# #     ax.scatter( X_bg, Z_bg, marker='o', s=5, c=DT_bg, cmap='RdBu_r', vmin=0, vmax=5e6 )

# ax.set_xlim(0,350)
# ax.set_ylim(5400,5650)

#
#
### PLOT PxRADIUS

lab_list = [
            '10 days',
            '23 days',
            '56 days',
            '4 months',
            '10 months',
            '23 months',
            '4 years',
            '10 years'
            ]


markers = ['o', 's', '^', 'v', 'd', 'p', '*', 'h']
cmap = plt.cm.jet
n_curves = len(selected_timesteps)
color_values = np.linspace(0.1, 0.92, 10)



# Select timetseps
t_min = 10*24*60*60
t_max = 10*365*24*60*60
log_times = np.logspace(np.log10(t_min), np.log10(t_max), 8)
selected_timesteps = []
for log_t in log_times:
    idx = np.argmin(np.abs(np.array(timesteps) - log_t))
    selected_timesteps.append(timesteps[idx])
selected_timesteps = sorted(set(selected_timesteps))
sel_strs = [ format_label(t) for t in selected_timesteps ] 
print(f"Plotting timesteps: {sel_strs}")

# Do the plotting
fig, ax = plt.subplots(figsize=(9/2.54, 6/2.54))
i=0
for t_bg in reversed(selected_timesteps) :
    lab = lab_list[-i-1]
    marker = markers[i % len(markers)]  # Cycle through markers
    color = cmap(color_values[-i])

    i+=1

    df_bg = df[df['t'] == t_bg]
    df_bg = df_bg[df['z'] == 5505]  # single depth
    
    X  = df_bg['x'].values
    P = df_bg['PRES'].values

    ax.plot(X, P/1e6, ls='--', lw=0.6, 
            label=lab, color=color)
#             marker=marker, markevery=(offset,1), markersize=3,markeredgewidth=0.2, fillstyle='none', 

#     ## ADD MARKERS
#     x_min, x_max = 0.2, 90
#     if i%2 : 
#         fac = 1.8
#         x_min *= fac
#         x_max *= fac
#     n_markers = 5  # Number of markers you want
#     x_marker_positions = np.geomspace(x_min, x_max, n_markers)
#     P_interpolated = np.interp(x_marker_positions, X, P)
#     ax.plot(x_marker_positions, P_interpolated/1e6, 
#             ls='', marker=marker, markersize=4,
#             markeredgewidth=0.3, #fillstyle='none', 
#             label=lab, markeredgecolor='k', color=color)    
#     print(f"x_markers:{x_marker_positions}  //  P_interpolated:{P_interpolated}")

ax.set_xlabel("Distance from the well (m)")
ax.set_ylabel(r"Pressure (MPa)")
ax.set_title(r"5m into the reservoir")
ax.set_xlim(0.1,100)
ax.set_xscale('log')
ax.set_ylim(63,75)
ax.legend( ncol=2, fontsize=6, loc='lower left' )

from matplotlib.ticker import FuncFormatter
ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f'{x:.10g}'))


fig.savefig('png/pres_profile.png', dpi=300)

plt.show()



