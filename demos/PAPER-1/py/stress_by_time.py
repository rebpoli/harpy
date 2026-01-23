#!/usr/bin/env -S python3 -i

import os
import subprocess
import shutil
import numpy as np
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count
import time
from timestr import format_time_duration
from netcdf import read_netcdf
from s3_helpers import (
    find_closest_xy_positions, extract_data_at_time, extract_vertical_profile,
    extract_time_series_at_position, find_closest_profile_for_depth, calculate_stress_ranges,
    create_profile_legend_elements, calculate_stress_differences, create_video_with_ffmpeg
)
from plot_util import savefig
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import FuncFormatter
def format_time(value, tick_number):
    if value < 1:
        if value < 0.1:
            return f''
        return f'{value:.1f}'  # 1 decimal place for values < 1
    else:
        return f'{int(value)}'  # No decimals for values >= 1



import matplotlib
matplotlib.use('Agg')

## Get the mplstyle
import os
cwd = os.path.dirname(__file__)
fn = f"{cwd}/my.mplstyle"
plt.style.use(fn)

# ==================== GLOBAL CONTROL VARIABLES ====================
NUM_LINES = 5  # Number of profile lines to display
TRACKED_DEPTHS = [5.0, -5.0]  # Depths to track (reservoir, caprock)
WELL_POSITION = (0.0, 0.0)  # Well position (x, y)
DEPTH_COLORS = ['red', 'blue']  # Colors for tracked depths
STRESS_COLORS = ['purple', 'orange']  # Colors for sigmaxx, sigmazz
XY_TOLERANCE = 1e-6  # Tolerance for (x,y) position matching
Z_TOLERANCE = 2.0  # Tolerance for depth matching (meters)

# Load dataset
filename = 'plane_yz.cd'
filepath = f"run/cdf/{filename}"
print(f"Loading {filepath}...")
dataset = read_netcdf(filepath)
print(dataset)

DIST_FROM_IFC = 5

def dataset_at( x,y,z ):
    x_coords = dataset.Coord.sel(vec3_comp='x')
    y_coords = dataset.Coord.sel(vec3_comp='y')
    z_coords = dataset.Coord.sel(vec3_comp='z')

    # Calculate Euclidean distance
    distances = np.sqrt((x_coords - x)**2 + (y_coords - y)**2 + (z_coords - z)**2)
    nearest_idx = distances.argmin().values
    nearest_point_idx = dataset.point_idx[nearest_idx].values

    print(f"Nearest point_idx: {nearest_point_idx}")
    print(f"Distance: {distances.min().values}")
    print(f"Coordinates: {dataset.Coord.isel(point_idx=nearest_idx).values}")

    return dataset.isel(point_idx=nearest_idx)

# Caprock
data_cap = dataset_at(0,0,-DIST_FROM_IFC)
df_cap = data_cap[["Total Stress","S3 Magnitude","Temperature"]].to_dataframe().reset_index()

# Resrevoir
data_res = dataset_at(0,0,DIST_FROM_IFC)
df_res = data_res[["Total Stress","S3 Magnitude","Pressure","Temperature"]].to_dataframe().reset_index()

# Save csv.
print("Save csv ...")
df_res.to_csv("csv/df_res.csv")
df_cap.to_csv("csv/df_cap.csv")
print("Done")

# import pandas as pd
# df_res = pd.read_csv("csv/df_res.csv")
# df_cap = pd.read_csv("csv/df_cap.csv")

## Time in days
df_res["time_d"] = np.sqrt(df_res.time/60/60/24)
df_cap["time_d"] = np.sqrt(df_cap.time/60/60/24)
df_res["time_d"] = df_res.time/60/60/24
df_cap["time_d"] = df_cap.time/60/60/24

# Do some plotting
fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(16/2.54,6/2.64), sharey=True)


## RES
dfxx = df_res[(df_res.ten9_comp=="xx") | (df_res.ten9_comp=="yy")]
s3 = df_res.groupby('time_d')[['time','S3 Magnitude']].max().reset_index()
sxx = dfxx.groupby('time_d')['Total Stress'].max().reset_index()
dfzz = df_res[df_res.ten9_comp=="zz"]
szz = dfzz.groupby('time_d')['Total Stress'].max().reset_index()
# s3_0 = s3[s3.time<0]['S3 Magnitude'].values[0]
s3_0 = df_res[df_res.time<0]['S3 Magnitude'].max()

ax1.axhline(-s3_0/1e6, ls='-', lw=1, alpha=1, c='orange', label=r'$\sigma_3(t=0)$')
ax1.plot(s3.time_d, -s3['S3 Magnitude']/1e6, c='k', alpha=0.3, lw=2, label=r"$\sigma_3$")
ax1.plot(sxx.time_d, -sxx['Total Stress']/1e6, c='blue',lw=0.7, ls='--', label=r"$\sigma_h$")
ax1.plot(szz.time_d, -szz['Total Stress']/1e6, c='green', lw=0.7, ls='-.', label=r"$\sigma_{zz}$")
ax1.plot(df_res.time_d, df_res.Pressure/1e6, c='red', lw=0.7, ls='-.', label=r"$p$")

## Caprock
dfxx = df_cap[df_cap.ten9_comp=="xx"]
s3 = df_cap.groupby('time_d')['S3 Magnitude'].max().reset_index()
sxx = dfxx.groupby('time_d')['Total Stress'].max().reset_index()
dfzz = df_cap[df_cap.ten9_comp=="zz"]
szz = dfzz.groupby('time_d')['Total Stress'].max().reset_index()
# s3_0 = s3[s3.time_d<0]['S3 Magnitude'].values[0]
s3_0 = df_cap[df_cap.time<0]['S3 Magnitude'].max()


ax2.axhline(-s3_0/1e6, ls='-', lw=1, alpha=1, c='orange', label=r'$\sigma_3^0$')
ax2.plot(s3.time_d, -s3['S3 Magnitude']/1e6, c='k', alpha=0.3, lw=2, label=r"$\sigma_3$")
ax2.plot(sxx.time_d, -sxx['Total Stress']/1e6, c='blue', lw=0.7, ls='--', label=r"$\sigma_h$")
ax2.plot(szz.time_d, -szz['Total Stress']/1e6, c='green', lw=0.7, ls='-.',label=r"$\sigma_{zz}$")

# Decorations
ax1.set_title(f"Near-well, {DIST_FROM_IFC}m into the reservoir")
ax2.set_title(f"Near-well, {DIST_FROM_IFC}m into the caprock")
for ax in [ax1,ax2] : 
#     ax2.set_xlabel(r"$\sqrt{\text{time}}$ (days)")
    ax.set_xlabel("Time (days)")
#     ax.set_xlim(0,np.sqrt(365*2))
#     ax.set_xlim(0,365*2)
    ax.set_xlim(0.05,365*3)
    ax.set_xscale('log')
    ax.xaxis.set_major_formatter(FuncFormatter(format_time))


lines1, labels1   = ax1.get_legend_handles_labels()
# leg=fig.legend(lines1, labels1, loc='lower center', bbox_to_anchor=(0.5,0.15), ncol=5,
#            fontsize=7, frameon=True, edgecolor='black', fancybox=False, framealpha=1)
leg=ax2.legend(lines1, labels1, loc='upper left', bbox_to_anchor=(1.05,1), ncol=1,
           fontsize=7, frameon=True, edgecolor='black', fancybox=False, framealpha=1)
leg.get_frame().set_linewidth(0.5)

ax1.set_ylabel("$-$Stress,  Pressure (MPa)")
# ax1.invert_yaxis()


fig.savefig("png/stress_by_time.png", dpi=500)
# plt.show()


#
#
# SHORT TERM PLOT
#
#

fig,[ ax1, ax2 ] = plt.subplots(2, 1, figsize=(9/2.54,7/2.64), sharex=True, gridspec_kw={'height_ratios': [2, 1]})
## RES
dfxx = df_res[(df_res.ten9_comp=="xx") | (df_res.ten9_comp=="yy")]
s3 = df_res.groupby('time_d')['S3 Magnitude'].max().reset_index()
sxx = dfxx.groupby('time_d')['Total Stress'].max().reset_index()
dfzz = df_res[df_res.ten9_comp=="zz"]
szz = dfzz.groupby('time_d')['Total Stress'].max().reset_index()
# s3_0 = s3[s3.time_d<0]['S3 Magnitude'].values[0]
s3_0 = df_res[df_res.time<0]['S3 Magnitude'].max()

ax1.plot(s3.time_d, -s3['S3 Magnitude']/1e6, c='k', alpha=0.5, lw=2, label=r"$\sigma_3$")
ax1.plot(sxx.time_d, -sxx['Total Stress']/1e6, c='blue',lw=0.7, ls='--', label=r"$\sigma_h$")
ax1.plot(szz.time_d, -szz['Total Stress']/1e6, c='green', lw=0.7, ls='-.', label=r"$\sigma_{zz}$")
ax1.plot(df_res.time_d, df_res.Pressure/1e6, c='red', lw=0.7, ls='-.', label=r"$p$")
ax1.axhline(-s3_0/1e6, ls='-', lw=0.9, alpha=0.8, c='orange', label=r'$\sigma_3(t=0)$')

ax2.plot(df_res.time_d, df_res.Temperature-273, c='#E0115F', lw=0.6, ls='--', label=r'$T$')
ax2.set_ylim(30,90)
ax2.set_ylabel(r"Temperature (°C)")

# Decorations
ax1.set_title(f"Near-well, {DIST_FROM_IFC}m into the reservoir")
# ax2.set_xlabel(r"$\sqrt{\text{time}}$ (days)")
ax2.set_xlabel(r"Time (days)")
ax1.set_xlim(0,np.sqrt(2))

ax1.set_ylabel(r"$-\sigma$, p (MPa)")
# ax1.invert_yaxis()

lines1, labels1   = ax1.get_legend_handles_labels()
leg = ax1.legend(lines1 , labels1 , 
                 loc='lower center', bbox_to_anchor=(0.55,0.02), ncol=6,
                 fontsize=7, frameon=True, edgecolor='black', fancybox=False,
                 framealpha=1, columnspacing=1)
leg.get_frame().set_linewidth(0.5)


lines1t, labels1t = ax2.get_legend_handles_labels()
leg = ax2.legend( lines1t,  labels1t, 
                 loc='upper center', bbox_to_anchor=(0.55,.98), ncol=6,
                 fontsize=7, frameon=True, edgecolor='black', fancybox=False,
                 framealpha=1, columnspacing=1)
leg.get_frame().set_linewidth(0.5)


print("Saving short term...")
fig.savefig("png/stress_by_time_short_term.png", dpi=500)
