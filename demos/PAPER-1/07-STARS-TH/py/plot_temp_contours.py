#!/usr/bin/env -S python3
"""
Plot DELTA_TEMP = 5 contours on XZ plots for selected timesteps.
DELTA_TEMP = TEMP - TEMP(t=0)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

from util import format_label

import os
cwd = os.path.dirname(__file__)
fn = f"{cwd}/my.mplstyle"
plt.style.use(fn)


# ============================================================================

# Option 4: Auto-select evenly spaced timesteps (default)
AUTO_SELECT = 12
contour_levels = [-5]

# ============================================================================

# Load the gzipped CSV file
print("Loading data...")
df = pd.read_csv('temperature.csv.gz', compression='gzip')

# Get initial temperature field at t=0
df_t0 = df[df['t'] == 0].copy()
df_t0 = df_t0.rename(columns={'TEMP': 'TEMP_0'})
df_t0 = df_t0[['x', 'z', 'TEMP_0']]

# We do not want very short noisy times
df = df[df.t > 0.1]

# Merge to get TEMP_0 for all timesteps
df = df.merge(df_t0, on=['x', 'z'], how='left')

# Calculate DELTA_TEMP
df['DELTA_TEMP'] = df['TEMP'] - df['TEMP_0']

# Get unique timesteps (excluding t=0 since DELTA_TEMP would be 0)
timesteps = sorted(df['t'].unique())
timesteps = [t for t in timesteps if t > 0]

print(f"\nTotal available timesteps: {len(timesteps)}")
print(f"Time range: {timesteps[0]:.2f} to {timesteps[-1]:.2e} days")

# Create logarithmically spaced time values
t_min = 10*24*60*60
t_max = 3*365*24*60*60
log_times = np.logspace(np.log10(t_min), np.log10(t_max), AUTO_SELECT)

selected_timesteps = []
for log_t in log_times:
    idx = np.argmin(np.abs(np.array(timesteps) - log_t))
    selected_timesteps.append(timesteps[idx])
selected_timesteps = sorted(set(selected_timesteps))

sel_strs = [ format_label(t) for t in selected_timesteps ] 

print(f"Plotting timesteps: {sel_strs}")
print(f"Contour levels: {contour_levels}")

# Create the plot
fig, ax1 = plt.subplots(1,1,figsize=(9/2.54, 6/2.54))

# Define colormap for different timesteps
colors = plt.cm.viridis(np.linspace(0, 1, len(selected_timesteps)))

lab_list = [
    [ '10 days', (10, 5550), 80 ],   
    [ ], 
    [ ], 
    [ ], 
    [ ], 
    [ ], 
    [ '4 months', (100, 5540), 55 ], 
    [ ], 
    [ ], 
    [ '16 months', (100, 5540), 38 ], 
    [ ], 
    [ '3 years', (150, 5550), 28 ], 
]

# Plot contours for each selected timestep
for idx, t in enumerate(selected_timesteps):
    df_t = df[df['t'] == t]
    
    # Get data
    x = df_t['x'].values
    z = df_t['z'].values
    delta_temp = df_t['DELTA_TEMP'].values
    
    print(f"\nTimestep {t:.2e}: DELTA_TEMP range = [{delta_temp.min():.2f}, {delta_temp.max():.2f}]")
    
    # Create a fine interpolation grid for smooth contours
    x_min, x_max = x.min(), x.max()
    z_min, z_max = z.min(), z.max()
    
    # Create a fine grid (increase resolution for smoother contours)
    n_points = 1000  # Number of points in each direction
    xi = np.logspace(np.log10(x_min), np.log10(x_max), n_points)
    zi = np.linspace(z_min, z_max, n_points)
    X, Z = np.meshgrid(xi, zi)
    
    # Interpolate DELTA_TEMP onto the fine grid
    DELTA_TEMP_grid = griddata((x, z), delta_temp, (X, Z), method='linear')
    
    # Check which contour levels are in range for this timestep
    valid_levels = [level for level in contour_levels 
                    if delta_temp.min() <= level <= delta_temp.max()]
    
    contour = ax1.contour(X, Z, DELTA_TEMP_grid, levels=valid_levels, 
                         colors='k', ls='--', linewidths=0.6)
        
    # Add inline labels to the contour
    l = lab_list[idx]
    if len(l) :
        labs = ax1.clabel(contour, inline=False, fontsize=6, 
                          manual=[l[1]], fmt=l[0], inline_spacing=0,
                          )#, 
        labs[0].set_rotation(l[2])
        x,y = labs[0].get_position()
        theta = np.deg2rad(l[2])
        nx=np.sin(theta)
        ny=np.cos(theta)
        offset=-5

        labs[0].set_position((x+offset*nx,y+offset*ny))


ax1.axhspan(5500, 5600, facecolor='gray', alpha=0.3, zorder=0)

# Labels and formatting
ax1.set_xlabel('Distance from well (m)')
ax1.set_ylabel('Depth (m)')
ax1.set_xlim(0,350)
ax1.set_ylim(5400,5650)
ax1.invert_yaxis()

# Create title based on contour levels
# if len(contour_levels) == 1:
#     title = f'DELTA_TEMP = {contour_levels[0]} Contours'
# else:
#     title = f'DELTA_TEMP Contours at levels: {contour_levels}'

ax1.set_title(r"Isothermal contours $\Delta\theta(t)=-5$°C, at selected times")

# Only show legend if we have contours
if ax1.get_legend_handles_labels()[0]:
    ax1.legend(loc='best')
ax1.grid(True)


fig.savefig('png/temperature_contours.png', dpi=300, bbox_inches='tight')



# Select timetseps
t_min = 10*24*60*60
t_max = 3*365*24*60*60
log_times = np.logspace(np.log10(t_min), np.log10(t_max), 7)
selected_timesteps = []
for log_t in log_times:
    idx = np.argmin(np.abs(np.array(timesteps) - log_t))
    selected_timesteps.append(timesteps[idx])
selected_timesteps = sorted(set(selected_timesteps))
sel_strs = [ format_label(t) for t in selected_timesteps ] 
print(f"Plotting timesteps: {sel_strs}")

lab_list = [
[ '10 days', 57, 87, 4], 
[], 
['45 days', 60, 82, 2 ], 
[], 
['8 months', 55, 75, 0], 
[], 
[ '3 years', 53, 57, -1 ], 
]


# Do the plotting
# fig, [ax1,ax2] = plt.subplots(1,2,figsize=(13/2.54, 6/2.54), sharey=True)
fig, ax1 = plt.subplots(figsize=(9/2.54, 6/2.54))
for idx, t_bg in enumerate(selected_timesteps) :
    df_bg = df[df['t'] == t_bg]
    df_bg = df_bg[df['z'] == 5505]  # single depth
    
    X  = df_bg['x'].values
    T = df_bg['TEMP'].values

    ax1.plot(X, T-273, ls='--', lw=0.6, color='k')

    ## Annotate
    l=lab_list[idx]
    if len(l) :
        trg_T = l[1] + 273 # target temperature
        txt = l[0] # annotation text
        theta = np.deg2rad(l[2]) # rotation angle
        print(f"Annotating '{txt}'")
        nx=-np.sin(theta)
        ny=np.cos(theta)
        offset=l[3]
        si = np.argmin(np.abs(T - trg_T))
        print(f"   si:{si}")
        x = X[si] + nx*offset
        y = T[si] -273 + ny*offset
        print(f"   (x,y):{x},{y}")
        ax1.text(x, y, txt, rotation=l[2], ha='center', va='bottom', fontsize=6)




ax1.set_xlabel("Distance from the well (m)")
ax1.set_ylabel(r"Temperature (°C)")
ax1.set_title(r"5m into the reservoir")
ax1.set_xlim(0,350)
# ax1.set_yscale('log')
# ax1.set_ylim(64,74)

fig.savefig('png/temp_profile_res.png', dpi=300)

#
#
# TEMP PROFILE CAPROCK
#

lab_list = [
[], 
[], 
[], 
['3 months', 69,  22, 0.3 ], 
['8 months', 65 , 22, 0.5], 
[ '17 months', 63, 25, 0.1], 
[ '3 years', 60,   25, 0 ], 
]

# Do the plotting
fig, ax2 = plt.subplots(figsize=(9/2.54, 6/2.54))
for idx, t_bg in enumerate(selected_timesteps) :
    df_bg = df[df['t'] == t_bg]
    df_bg = df_bg[df['z'] == 5495]  # single depth
    
    X  = df_bg['x'].values
    T = df_bg['TEMP'].values

    ax2.plot(X, T-273, ls='--', lw=0.6, color='k')

    ## Annotate
    l=lab_list[idx]
    if len(l) :
        trg_T = l[1] + 273 # target temperature
        txt = l[0] # annotation text
        theta = np.deg2rad(l[2]) # rotation angle
        print(f"Annotating '{txt}'")
        nx=-np.sin(theta)
        ny=np.cos(theta)
        offset=l[3]
        si = np.argmin(np.abs(T - trg_T))
        print(f"   si:{si}")
        x = X[si] + nx*offset
        y = T[si] -273 + ny*offset
        print(f"   (x,y):{x},{y}")
        ax2.text(x, y, txt, rotation=l[2], ha='center', va='bottom', fontsize=6, rotation_mode='anchor')

ax2.set_xlabel("Distance from the well (m)")
ax2.set_ylabel(r"Temperature (°C)")
ax2.set_xlim(0,350)
ax2.set_title(r"5m into the caprock")
# ax2.set_yscale('log')
# ax2.set_ylim(64,74)

fig.savefig('png/temp_profile_caprock.png', dpi=300)



plt.show()
