#!/usr/bin/env -S python 

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.interpolate import griddata
import time
from plot_util import savefig

import matplotlib
matplotlib.use('Agg')

from netcdf import read_netcdf
from timestr import format_time_duration

#
#
def mohr_coulomb_envelope(phi, cohesion, pmax, npoints=200):
    p = np.linspace(0, pmax, npoints)
    q = p * np.tan(phi) + cohesion / np.cos(phi)
    C1 = 6 * cohesion * np.cos(phi) / ( 3 - np.sin(phi)
    C2 = 6 * np.sin(phi) / ( 3 - np.sin(phi)
    q = C1 + C2 * p
    return -p, q

## Get the mplstyle
import os
cwd = os.path.dirname(__file__)
fn = f"{cwd}/my.mplstyle"
plt.style.use(fn)

# ---------------------------------------------------------------------
#          CONFIGURABLE PARAMETERS - PLASTICITY ENVELOPE
# ---------------------------------------------------------------------
phi_salt_deg = 40.0       # Friction angle for salt [degrees]
phi_carb_deg = 50.0       # Friction angle for carbonate [degrees]
cohesion_salt = 2.0       # Cohesion for salt [MPa]
cohesion_carb = 10.0      # Cohesion for carbonate [MPa]
p_envelope_max = 80.0     # Maximum p' to draw the envelope up to [MPa]
phi_salt = np.radians(phi_salt_deg)
phi_carb = np.radians(phi_carb_deg)

# READ DATASET AND GRAB SOME HANDLERS
filepath = f"run/cdf/plane_yz.cd"
ds = read_netcdf(filepath)
ds['Delta S3'] = ds['S3 Magnitude'] - ds['S3 Magnitude'].isel(time=0)
time_values = ds.time.values
#   Coords setup
coords = ds['Coord'].values
x_comp_idx = 1 # X axis=Y (1)
y_comp_idx = 2 # Y axis=Z (2)
vector_density = 15
coords_2d = coords[:, [x_comp_idx, y_comp_idx]]
x_min, x_max = coords_2d[:, 0].min(), coords_2d[:, 0].max()
y_min, y_max = coords_2d[:, 1].min(), coords_2d[:, 1].max()

# Extract the variables
inv_p = ds['Invariant P(eff)'].values
inv_q = ds['Invariant Q'].values
time = ds['time'].values
z_coord = ds['Coord'].values[:, 2]  # Z is the 3rd component (index 2)

# Flatten the arrays
inv_p_flat = inv_p.flatten() / 1e6
inv_q_flat = inv_q.flatten() / 1e6
n_points = inv_p.shape[1]
time_flat = np.repeat(time, n_points)
z_flat = np.tile(z_coord, len(time))  # Repeat Z coords for each time step

# Convert time to days
time_flat_days = time_flat / 86400.0

# Subsample the data
sample_fraction = 1
n_samples = int(len(inv_p_flat) * sample_fraction)
indices = np.random.choice(len(inv_p_flat), size=n_samples, replace=False)

inv_p_sub = inv_p_flat[indices]
inv_q_sub = inv_q_flat[indices]
time_sub = time_flat_days[indices]
z_sub = z_flat[indices]

# Split into reservoir (Z>0) and caprock (Z<0)
reservoir_mask = z_sub > 0
caprock_mask = z_sub < 0

# Convert 15cm x 8cm to inches
fig = plt.figure(figsize=(15/2.54, 6.5/2.54))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.08], wspace=0.02, left=0.075, right=0.92,top=0.92,bottom=0.12)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax_cbar = fig.add_subplot(gs[0, 2])

alpha=0.1
cmap = "jet"

## Colorbar
vmin=0
vmax = time_flat_days.max()
import matplotlib.colors as mcolors
norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=(vmin+vmax)/2, vmax=vmax)

# Plot caprock (Z<0)
scatter1 = ax1.scatter(inv_p_sub[caprock_mask], inv_q_sub[caprock_mask],
                       c=time_sub[caprock_mask],
                       cmap=cmap,
                       alpha=alpha,
                       norm=norm,
                       s=3,
                       edgecolors='none')

p_env, q_env = mohr_coulomb_envelope(phi_salt, cohesion_salt, p_envelope_max)
ax1.plot(p_env, q_env, color='k', linestyle="--", linewidth=0.4, label="MC envelope")
ax1.fill_between( p_env, q_env, q_env.max(), alpha=0.2, color='k') 

# Plot reservoir (Z>0)
scatter2 = ax2.scatter(inv_p_sub[reservoir_mask], inv_q_sub[reservoir_mask],
                       c=time_sub[reservoir_mask],
                       cmap=cmap,
                       norm=norm,
                       alpha=alpha,
                       s=3,
                       edgecolors='none')

p_env, q_env= mohr_coulomb_envelope(phi_carb, cohesion_carb, p_envelope_max)
ax2.plot(p_env, q_env, color='k', linestyle="--", linewidth=0.4, label="MC envelope")
ax2.fill_between( p_env, q_env, q_env.max(), alpha=0.2, color='k') 

cbar = plt.colorbar(scatter1, cax=ax_cbar, drawedges=False)
cbar.solids.set_alpha(1.0)
cbar.set_label(r'Time (days)', rotation=270, labelpad=10)

ax1.set_title('Caprock')
ax2.set_title('Reservoir')

for ax in [ ax1, ax2 ] :
    ax.set_xlabel("$p'$ (MPa)")
    ax.set_ylabel("$q$ (MPa)")
    ax.set_xlim(-40,0)
    ax.set_ylim(0,30)
    ax.invert_xaxis()

savefig( fig, "png/PQ_contact.png" )

