#!/usr/bin/env -S python -i

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
from custom_cmap import create_custom_cmap, create_custom_cmap_nozero, create_custom_cmap_nomin

# Set global font properties
plt.rcParams['font.family'] = 'Cambria Math'
plt.rcParams['font.size'] = 7
plt.rcParams['axes.titlesize'] = 8
plt.rcParams['axes.labelsize'] = 7
plt.rcParams['xtick.labelsize'] = 7
plt.rcParams['ytick.labelsize'] = 7
plt.rcParams['legend.fontsize'] = 7
plt.rcParams['figure.titlesize'] = 9
plt.rcParams['axes.labelpad'] = 1
plt.rcParams['axes.titlepad'] = 4
plt.rcParams['legend.edgecolor'] = 'none'
plt.rcParams['legend.framealpha'] = 1
plt.rcParams['legend.handlelength'] = 0.8



# plt.rcParams['lines.linewidth'] = 2
# plt.rcParams['axes.linewidth'] = 1.5
# plt.rcParams['grid.alpha'] = 0.3



#
#
#
def create_vector_segments(coords, vectors, pt_length=10):
    (center_x, center_y), (vec_x, vec_y) = coords[:, :2].T, vectors[:, :2].T
    return np.column_stack([center_x, center_y, vec_x, vec_y, np.full(len(center_x), pt_length)])

#
#
#
def subsample_data_indices(coords_2d, target_density):
    x_coords, y_coords = coords_2d[:, 0], coords_2d[:, 1]
    x_min, x_max = x_coords.min(), x_coords.max()
    y_min, y_max = y_coords.min(), y_coords.max()
    
    x_bins = np.linspace(x_min, x_max, target_density + 1)
    y_bins = np.linspace(y_min, y_max, target_density + 1)
    
    selected_indices = []
    for i in range(target_density):
        for j in range(target_density):
            x_mask = (x_coords >= x_bins[i]) & (x_coords < x_bins[i+1])
            y_mask = (y_coords >= y_bins[j]) & (y_coords < y_bins[j+1])
            candidates = np.where(x_mask & y_mask)[0]
            if len(candidates) > 0:
                selected_indices.append(candidates[0])
    
    return np.array(selected_indices)

#
#
#
def plot_segments(ax, segments, color, lengthscale, linewidth, alpha):
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    fig = ax.get_figure()
    fig_width_inch, fig_height_inch = fig.get_size_inches()
    
    bbox = ax.get_position()
    axis_width_points = bbox.width * fig_width_inch * 72
    axis_height_points = bbox.height * fig_height_inch * 72
    
    data_width, data_height = xlim[1] - xlim[0], ylim[1] - ylim[0]
    scale_x = data_width / axis_width_points
    scale_y = data_height / axis_height_points
    
    for center_x, center_y, vec_x, vec_y, pt_length in segments:
        if not (np.isnan(center_x) or np.isnan(center_y) or np.isnan(vec_x) or np.isnan(vec_y)):
            data_length_x = pt_length * scale_x
            data_length_y = pt_length * scale_y
            
            dx = vec_x * data_length_x / 2 * lengthscale
            dy = vec_y * data_length_y / 2 * lengthscale
            
            x_coords = [center_x - dx, center_x + dx]
            y_coords = [center_y - dy, center_y + dy]
            
            line = Line2D(x_coords, y_coords, color=color, linewidth=linewidth, alpha=alpha)
            ax.add_line(line)

def plot_s1_s3(ax, s1_seg, s3_seg) :
    plot_segments(ax, s3_seg, 'lightgreen', 0.35, 0.7, 1)
    plot_segments(ax, s1_seg, 'black',      0.6, 0.7, 0.8)

#
#
#
#
#


# READ DATASET AND GRAB SOME HANDLERS
filepath = f"run/cdf/plane_yz.cd"
ds = read_netcdf(filepath)
ds['Delta S3'] = ds['S3 Magnitude'] - ds['S3 Magnitude'][0,:]

print(ds)
time_values = ds.time.values
#   Coords setup
coords = ds['Coord'].values
x_comp_idx = 1 # X axis=Y (1)
y_comp_idx = 2 # Y axis=Z (2)
vector_density = 15
coords_2d = coords[:, [x_comp_idx, y_comp_idx]]
x_min, x_max = coords_2d[:, 0].min(), coords_2d[:, 0].max()
y_min, y_max = coords_2d[:, 1].min(), coords_2d[:, 1].max()

#
# PLOT - VP STRAIN 
#

# CREATE SUBPLOT AND AXES
fig = plt.figure(figsize=(15/2.54, 6.5/2.54))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.08], wspace=0.1, left=0.075, right=0.92,top=0.92,bottom=0.12)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax_cbar = fig.add_subplot(gs[0, 2])

# SETUP COLORBAR
# vmin = -1e-12
# vmax = 1e-12
_data = -ds['VP Strain'].sel(ten9_comp='zz')
vmin, vmax = float(_data.quantile(0.05).values), float(_data.quantile(0.95).values)
vzero = 0 #0.5 * (vmin+vmax)
cmap = create_custom_cmap(vmin,vzero,vmax)
dummy_data = np.linspace(vmin, vmax, 100).reshape(10, 10) * 1e4
im = ax1.imshow(dummy_data, cmap=cmap, vmin=vmin*1e4, vmax=vmax*1e4)
ax1.clear()
cbar = plt.colorbar(im, cax=ax_cbar)
cbar.set_label(r'Viscoplastic vertical strain ($-\varepsilon_{\text{v}}^{\text{vp}} \text{}{}_{\times 10^{4}}$)', rotation=270, labelpad=10)

ax_cbar.annotate('', xy=(0.7, 0.5), xytext=(0.7, 0.3), xycoords='axes fraction', arrowprops=dict(arrowstyle='->', lw=0.5, color='black'))
ax_cbar.text(0.3, 0.4, 'Compressed', transform=ax_cbar.transAxes, rotation=90, va='center', ha='center', fontsize=6)
ax_cbar.annotate('', xy=(0.7, 0.05), xytext=(0.7, 0.25), xycoords='axes fraction', arrowprops=dict(arrowstyle='->', lw=0.5, color='black'))
ax_cbar.text(0.3, 0.15, 'Elongated', transform=ax_cbar.transAxes, rotation=90, va='center', ha='center', fontsize=6)

# PLOT
for time_idx, ax in zip( [32, 113 ] , [ax1, ax2 ] ) :
    print(f"Processing time {format_time_duration(time_values[time_idx])}")
    subsample_indices = subsample_data_indices(coords_2d, vector_density)
    s1_vectors = ds['S1'][time_idx, :, :].values[subsample_indices]
    s3_vectors = ds['S3'][time_idx, :, :].values[subsample_indices]
    s1_2d = s1_vectors[:, [x_comp_idx, y_comp_idx]]
    s3_2d = s3_vectors[:, [x_comp_idx, y_comp_idx]]
    subsampled_coords_2d = coords[subsample_indices][:, [x_comp_idx, y_comp_idx]]
    s1_segments = create_vector_segments(subsampled_coords_2d, s1_2d)
    s3_segments = create_vector_segments(subsampled_coords_2d, s3_2d)

    background_data = -ds['VP Strain'].sel(ten9_comp='zz')[time_idx, :].values
    bg_xi = np.linspace(x_min, x_max, 100)
    bg_yi = np.linspace(y_min, y_max, 100)
    all_coords_2d = coords[:, [x_comp_idx, y_comp_idx]]
    bg_Xi, bg_Yi = np.meshgrid(bg_xi, bg_yi)
    bg_Zi = griddata(all_coords_2d, background_data, (bg_Xi, bg_Yi), method='linear')

    # Do the plotting
    ax.imshow(bg_Zi, extent=[bg_Xi.min(), bg_Xi.max(), bg_Yi.min(), bg_Yi.max()],
             origin='lower', aspect='auto', interpolation='bilinear',
             cmap=cmap, vmin=vmin, vmax=vmax, alpha=0.6)
    plot_s1_s3( ax, s1_segments, s3_segments )

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.invert_yaxis()
    legend_elements = [ Line2D([0], [0], color='black', lw=1, label=r'$\sigma_1$'),
                        Line2D([0], [0], color='lightgreen', lw=1, label=r'$\sigma_3$') ]
    ax.legend(handles=legend_elements, loc='upper right')
    ax.set_xlabel("Distance from the well (m)")


ax2.tick_params(axis='y', labelleft=False)
ax1.set_ylabel("Distance from the caprock interface (m)")
ax1.set_title("1 month of injection")
ax2.set_title("1 year of injection")

savefig(fig, "png/PAPER-STRAIN-YZ.png" )

#
# PLOT - VP STRAIN RATE
#

# CREATE SUBPLOT AND AXES
fig = plt.figure(figsize=(15/2.54, 6.5/2.54))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.08], wspace=0.1, left=0.075, right=0.92,top=0.92,bottom=0.12)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax_cbar = fig.add_subplot(gs[0, 2])

# SETUP COLORBAR
# vmin = -1e-12
# vmax = 1e-12
_data = -ds['VP Strain Rate'].sel(ten9_comp='zz')
vmin, vmax = float(_data.quantile(0.05).values), float(_data.quantile(0.95).values)
vzero = 0 #0.5 * (vmin+vmax)
cmap = create_custom_cmap(vmin,vzero,vmax)
dummy_data = np.linspace(vmin, vmax, 100).reshape(10, 10)
im = ax1.imshow(dummy_data, cmap=cmap, vmin=vmin, vmax=vmax)
ax1.clear()
cbar = plt.colorbar(im, cax=ax_cbar)
cbar.set_label(r'Viscoplastic vertical strain rate ($-\dot{\varepsilon}_{v}^{\text{vp}}$)', rotation=270, labelpad=10)

ax_cbar.annotate('', xy=(0.7, 0.55), xytext=(0.7, 0.35), xycoords='axes fraction', arrowprops=dict(arrowstyle='->', lw=0.5, color='black'))
ax_cbar.text(0.3, 0.45, 'Compressing', transform=ax_cbar.transAxes, rotation=90, va='center', ha='center', fontsize=6)
ax_cbar.annotate('', xy=(0.7, 0.08), xytext=(0.7, 0.28), xycoords='axes fraction', arrowprops=dict(arrowstyle='->', lw=0.5, color='black'))
ax_cbar.text(0.3, 0.18, 'Elongating', transform=ax_cbar.transAxes, rotation=90, va='center', ha='center', fontsize=6)


# PLOT
for time_idx, ax in zip( [32, 113 ] , [ax1, ax2 ] ) :
    print(f"Processing time {format_time_duration(time_values[time_idx])}")
    subsample_indices = subsample_data_indices(coords_2d, vector_density)
    s1_vectors = ds['S1'][time_idx, :, :].values[subsample_indices]
    s3_vectors = ds['S3'][time_idx, :, :].values[subsample_indices]
    s1_2d = s1_vectors[:, [x_comp_idx, y_comp_idx]]
    s3_2d = s3_vectors[:, [x_comp_idx, y_comp_idx]]
    subsampled_coords_2d = coords[subsample_indices][:, [x_comp_idx, y_comp_idx]]
    s1_segments = create_vector_segments(subsampled_coords_2d, s1_2d)
    s3_segments = create_vector_segments(subsampled_coords_2d, s3_2d)

    background_data = -ds['VP Strain Rate'].sel(ten9_comp='zz')[time_idx,:].values
    bg_xi = np.linspace(x_min, x_max, 100)
    bg_yi = np.linspace(y_min, y_max, 100)
    all_coords_2d = coords[:, [x_comp_idx, y_comp_idx]]
    bg_Xi, bg_Yi = np.meshgrid(bg_xi, bg_yi)
    bg_Zi = griddata(all_coords_2d, background_data, (bg_Xi, bg_Yi), method='linear')

    # Do the plotting
    ax.imshow(bg_Zi, extent=[bg_Xi.min(), bg_Xi.max(), bg_Yi.min(), bg_Yi.max()],
             origin='lower', aspect='auto', interpolation='bilinear',
             cmap=cmap, vmin=vmin, vmax=vmax, alpha=0.6)
    plot_s1_s3( ax, s1_segments, s3_segments )

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.invert_yaxis()
    legend_elements = [ Line2D([0], [0], color='black', lw=1, label=r'$\sigma_1$'),
                        Line2D([0], [0], color='lightgreen', lw=1, label=r'$\sigma_3$') ]
    ax.legend(handles=legend_elements, loc='upper right')
    ax.set_xlabel("Distance from the well (m)")


ax2.tick_params(axis='y', labelleft=False)
ax1.set_ylabel("Distance from the caprock interface (m)")
ax1.set_title("1 month of injection")
ax2.set_title("1 year of injection")

savefig(fig, "png/PAPER-STRAIN-RATE-YZ.png" )

#
# PLOT - DELTA S3
#

# CREATE SUBPLOT AND AXES
fig = plt.figure(figsize=(15/2.54, 6.5/2.54))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.08], wspace=0.1, left=0.08, right=0.92,top=0.92,bottom=0.14)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax_cbar = fig.add_subplot(gs[0, 2])

# SETUP COLORBAR
vmin = -8
vmax = 2
dummy_data = np.linspace(vmin, vmax, 100).reshape(10, 10)
import matplotlib.colors as mcolors
cmap = create_custom_cmap(vmin,0,vmax)
im = ax1.imshow(dummy_data, cmap=cmap, vmin=vmin, vmax=vmax)
cbar = plt.colorbar(im, cax=ax_cbar)
cbar.set_label(r'Min. principal stress change ($-\Delta \sigma_3$, MPa)', rotation=270, labelpad=7, fontsize=7)
ax1.clear()

ax_cbar.annotate('', xy=(0.5, 0.85), xytext=(0.5, 0.67), xycoords='axes fraction', arrowprops=dict(arrowstyle='->', lw=0.5, color='black'))
ax_cbar.annotate('', xy=(0.5, 0.32), xytext=(0.5, 0.15),  xycoords='axes fraction', arrowprops=dict(arrowstyle='-', lw=0.5, color='black'))
ax_cbar.text(0.5, 0.5, 'More compressive', transform=ax_cbar.transAxes, rotation=90, va='center', ha='center', fontsize=6)

# PLOT
for time_idx, ax in zip( [32, 113 ] , [ax1, ax2 ] ) :
    print(f"Processing time {format_time_duration(time_values[time_idx])}")
    subsample_indices = subsample_data_indices(coords_2d, vector_density)
    s1_vectors = ds['S1'][time_idx, :, :].values[subsample_indices]
    s3_vectors = ds['S3'][time_idx, :, :].values[subsample_indices]
    s1_2d = s1_vectors[:, [x_comp_idx, y_comp_idx]]
    s3_2d = s3_vectors[:, [x_comp_idx, y_comp_idx]]
    subsampled_coords_2d = coords[subsample_indices][:, [x_comp_idx, y_comp_idx]]
    s1_segments = create_vector_segments(subsampled_coords_2d, s1_2d)
    s3_segments = create_vector_segments(subsampled_coords_2d, s3_2d)

    background_data = -ds['Delta S3'][time_idx, :].values / 1e6
    bg_xi = np.linspace(x_min, x_max, 100)
    bg_yi = np.linspace(y_min, y_max, 100)
    all_coords_2d = coords[:, [x_comp_idx, y_comp_idx]]
    bg_Xi, bg_Yi = np.meshgrid(bg_xi, bg_yi)
    bg_Zi = griddata(all_coords_2d, background_data, (bg_Xi, bg_Yi), method='linear')

    # Do the plotting
    ax.imshow(bg_Zi, extent=[bg_Xi.min(), bg_Xi.max(), bg_Yi.min(), bg_Yi.max()],
             origin='lower', aspect='auto', interpolation='bilinear',
             cmap=cmap, vmin=vmin, vmax=vmax, alpha=0.6)
    plot_s1_s3( ax, s1_segments, s3_segments)

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.invert_yaxis()
    legend_elements = [ Line2D([0], [0], color='black', lw=1, label=r'$\sigma_1$'),
                        Line2D([0], [0], color='lightgreen', lw=1, label=r'$\sigma_3$') ]
    ax.legend(handles=legend_elements, loc='upper right')
    ax.set_xlabel("Distance from the well (m)")


ax2.tick_params(axis='y', labelleft=False)
ax1.set_ylabel("Distance from the caprock interface (m)")
ax1.set_title("1 month of injection")
ax2.set_title("1 year of injection")

savefig(fig, "png/PAPER-S3-YZ.png" )

#
# PLOT - TEMPERATURE 
#

# CREATE SUBPLOT AND AXES
fig = plt.figure(figsize=(15/2.54, 6.5/2.54))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.08], wspace=0.1, left=0.08, right=0.92,top=0.92,bottom=0.14)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax_cbar = fig.add_subplot(gs[0, 2])

# SETUP COLORBAR
# vmin = -2
# vmax = 8
_data = ds['Temperature'] - 273
vmin, vmax = float(_data.quantile(0.05).values), float(_data.quantile(0.95).values)
vzero = (vmin+vmax)/2
dummy_data = np.linspace(vmin, vmax, 100).reshape(10, 10)
import matplotlib.colors as mcolors
cmap = create_custom_cmap(vmin,vzero,vmax)
im = ax1.imshow(dummy_data, cmap=cmap, vmin=vmin, vmax=vmax)
cbar = plt.colorbar(im, cax=ax_cbar)
cbar.set_label(r'Temperature (°C)', rotation=270, labelpad=7)
ax1.clear()

# PLOT
for time_idx, ax in zip( [32, 113 ] , [ax1, ax2 ] ) :
    print(f"Processing time {format_time_duration(time_values[time_idx])}")
    subsample_indices = subsample_data_indices(coords_2d, vector_density)
    s1_vectors = ds['S1'][time_idx, :, :].values[subsample_indices]
    s3_vectors = ds['S3'][time_idx, :, :].values[subsample_indices]
    s1_2d = s1_vectors[:, [x_comp_idx, y_comp_idx]]
    s3_2d = s3_vectors[:, [x_comp_idx, y_comp_idx]]
    subsampled_coords_2d = coords[subsample_indices][:, [x_comp_idx, y_comp_idx]]
    s1_segments = create_vector_segments(subsampled_coords_2d, s1_2d)
    s3_segments = create_vector_segments(subsampled_coords_2d, s3_2d)

    background_data = ds['Temperature'][time_idx, :].values -273
    bg_xi = np.linspace(x_min, x_max, 100)
    bg_yi = np.linspace(y_min, y_max, 100)
    all_coords_2d = coords[:, [x_comp_idx, y_comp_idx]]
    bg_Xi, bg_Yi = np.meshgrid(bg_xi, bg_yi)
    bg_Zi = griddata(all_coords_2d, background_data, (bg_Xi, bg_Yi), method='linear')

    # Do the plotting
    ax.imshow(bg_Zi, extent=[bg_Xi.min(), bg_Xi.max(), bg_Yi.min(), bg_Yi.max()],
             origin='lower', aspect='auto', interpolation='bilinear',
             cmap=cmap, vmin=vmin, vmax=vmax, alpha=0.6)
    plot_s1_s3( ax, s1_segments, s3_segments)

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.invert_yaxis()
    legend_elements = [ Line2D([0], [0], color='black', lw=1, label=r'$\sigma_1$'),
                        Line2D([0], [0], color='lightgreen', lw=1, label=r'$\sigma_3$') ]
    ax.legend(handles=legend_elements, loc='upper right')
    ax.set_xlabel("Distance from the well (m)")


ax2.tick_params(axis='y', labelleft=False)
ax1.set_ylabel("Distance from the caprock interface (m)")
ax1.set_title("1 month of injection")
ax2.set_title("1 year of injection")

savefig(fig, "png/PAPER-TEMP-YZ.png" )


#
# PLOT - PRESSURE
#

# CREATE SUBPLOT AND AXES
fig = plt.figure(figsize=(15/2.54, 6.5/2.54))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.08], wspace=0.1, left=0.08, right=0.92,top=0.92,bottom=0.14)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax_cbar = fig.add_subplot(gs[0, 2])

# SETUP COLORBAR
# vmin = -2
# vmax = 8
_data = ds['Pressure']/1e6
vmin, vmax = float(_data.quantile(0.05).values), float(_data.quantile(0.95).values)
vmin, vmax = ( 60, 67 )
# vzero = (vmin+vmax)/2
dummy_data = np.linspace(vmin, vmax, 100).reshape(10, 10)
import matplotlib.colors as mcolors
cmap = create_custom_cmap_nozero()
im = ax1.imshow(dummy_data, cmap=cmap, vmin=vmin, vmax=vmax)
cbar = plt.colorbar(im, cax=ax_cbar)
cbar.set_label(r'Pressure (MPa)', rotation=270, labelpad=7)
ax1.clear()

# PLOT
for time_idx, ax in zip( [32, 113 ] , [ax1, ax2 ] ) :
    print(f"Processing time {format_time_duration(time_values[time_idx])}")
    subsample_indices = subsample_data_indices(coords_2d, vector_density)
    s1_vectors = ds['S1'][time_idx, :, :].values[subsample_indices]
    s3_vectors = ds['S3'][time_idx, :, :].values[subsample_indices]
    s1_2d = s1_vectors[:, [x_comp_idx, y_comp_idx]]
    s3_2d = s3_vectors[:, [x_comp_idx, y_comp_idx]]
    subsampled_coords_2d = coords[subsample_indices][:, [x_comp_idx, y_comp_idx]]
    s1_segments = create_vector_segments(subsampled_coords_2d, s1_2d)
    s3_segments = create_vector_segments(subsampled_coords_2d, s3_2d)

    background_data = ds['Pressure'][time_idx, :].values/1e6
    bg_xi = np.linspace(x_min, x_max, 100)
    bg_yi = np.linspace(y_min, y_max, 100)
    all_coords_2d = coords[:, [x_comp_idx, y_comp_idx]]
    bg_Xi, bg_Yi = np.meshgrid(bg_xi, bg_yi)
    bg_Zi = griddata(all_coords_2d, background_data, (bg_Xi, bg_Yi), method='linear')

    # Do the plotting
    ax.imshow(bg_Zi, extent=[bg_Xi.min(), bg_Xi.max(), bg_Yi.min(), bg_Yi.max()],
             origin='lower', aspect='auto', interpolation='bilinear',
             cmap=cmap, vmin=vmin, vmax=vmax, alpha=0.6)
    plot_s1_s3( ax, s1_segments, s3_segments)

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.invert_yaxis()
    legend_elements = [ Line2D([0], [0], color='black', lw=1, label=r'$\sigma_1$'),
                        Line2D([0], [0], color='lightgreen', lw=1, label=r'$\sigma_3$') ]
    ax.legend(handles=legend_elements, loc='upper right')
    ax.set_xlabel("Distance from the well (m)")


ax2.tick_params(axis='y', labelleft=False)
ax1.set_ylabel("Distance from the caprock interface (m)")
ax1.set_title("1 month of injection")
ax2.set_title("1 year of injection")

savefig(fig, "png/PAPER-PRES-YZ.png" )

#
# PLOT - PRESSURE
#

# CREATE SUBPLOT AND AXES
fig = plt.figure(figsize=(15/2.54, 6.5/2.54))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.08], wspace=0.1, left=0.08, right=0.92,top=0.92,bottom=0.14)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax_cbar = fig.add_subplot(gs[0, 2])

# SETUP COLORBAR
# vmin = -2
# vmax = 8
_data = ds['Invariant Q']/1e6
vmin, vmax = float(_data.quantile(0.05).values), float(_data.quantile(0.95).values)
vmin, vmax = ( 0, 15 )
dummy_data = np.linspace(vmin, vmax, 100).reshape(10, 10)
import matplotlib.colors as mcolors
cmap = create_custom_cmap_nomin()
im = ax1.imshow(dummy_data, cmap=cmap, vmin=vmin, vmax=vmax)
cbar = plt.colorbar(im, cax=ax_cbar)
cbar.set_label(r"Invariant Q' (MPa)", rotation=270, labelpad=7)
ax1.clear()

# PLOT
for time_idx, ax in zip( [32, 113 ] , [ax1, ax2 ] ) :
    print(f"Processing time {format_time_duration(time_values[time_idx])}")
    subsample_indices = subsample_data_indices(coords_2d, vector_density)
    s1_vectors = ds['S1'][time_idx, :, :].values[subsample_indices]
    s3_vectors = ds['S3'][time_idx, :, :].values[subsample_indices]
    s1_2d = s1_vectors[:, [x_comp_idx, y_comp_idx]]
    s3_2d = s3_vectors[:, [x_comp_idx, y_comp_idx]]
    subsampled_coords_2d = coords[subsample_indices][:, [x_comp_idx, y_comp_idx]]
    s1_segments = create_vector_segments(subsampled_coords_2d, s1_2d)
    s3_segments = create_vector_segments(subsampled_coords_2d, s3_2d)

    background_data = ds['Invariant Q'][time_idx, :].values/1e6
    bg_xi = np.linspace(x_min, x_max, 100)
    bg_yi = np.linspace(y_min, y_max, 100)
    all_coords_2d = coords[:, [x_comp_idx, y_comp_idx]]
    bg_Xi, bg_Yi = np.meshgrid(bg_xi, bg_yi)
    bg_Zi = griddata(all_coords_2d, background_data, (bg_Xi, bg_Yi), method='linear')

    # Do the plotting
    ax.imshow(bg_Zi, extent=[bg_Xi.min(), bg_Xi.max(), bg_Yi.min(), bg_Yi.max()],
             origin='lower', aspect='auto', interpolation='bilinear',
             cmap=cmap, vmin=vmin, vmax=vmax, alpha=0.6)
    plot_s1_s3( ax, s1_segments, s3_segments)

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.invert_yaxis()
    legend_elements = [ Line2D([0], [0], color='black', lw=1, label=r'$\sigma_1$'),
                        Line2D([0], [0], color='lightgreen', lw=1, label=r'$\sigma_3$') ]
    ax.legend(handles=legend_elements, loc='upper right')
    ax.set_xlabel("Distance from the well (m)")


ax2.tick_params(axis='y', labelleft=False)
ax1.set_ylabel("Distance from the caprock interface (m)")
ax1.set_title("1 month of injection")
ax2.set_title("1 year of injection")

savefig(fig, "png/PAPER-Q-INV-YZ.png" )
