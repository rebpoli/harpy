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


z = dataset['Coord'].sel(vec3_comp= 'z').values

# Caprock
data_cap = dataset.isel(point_idx = z==-5)
df_cap = data_cap[["Total Stress","S3 Magnitude"]].to_dataframe().reset_index()

# Resrevoir
data_res = dataset.isel(point_idx = z==5)
df_res = data_res[["Total Stress","S3 Magnitude"]].to_dataframe().reset_index()

# Save csv.
df_res.to_csv("csv/df_res.csv")
df_cap.to_csv("csv/df_cap.csv")

import pandas as pd
df_res = pd.read_csv("csv/df_res.csv")
df_cap = pd.read_csv("csv/df_cap.csv")

## Time in days
df_res["time_d"] = df_res.time/60/60/24
df_cap["time_d"] = df_cap.time/60/60/24

# Do some plotting
fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(17.5/2.54,6/2.64), sharey=True)


## RES
dfxx = df_res[df_res.ten9_comp=="xx"]
s3 = df_res.groupby('time_d')['S3 Magnitude'].max().reset_index()
sxx = dfxx.groupby('time_d')['Total Stress'].max().reset_index()
dfzz = df_res[df_res.ten9_comp=="zz"]
szz = dfzz.groupby('time_d')['Total Stress'].max().reset_index()
s3_0 = s3[s3.time_d<0]['S3 Magnitude'].values[0]

ax1.axhline(s3_0/1e6, ls=':', lw=2, alpha=0.3, c='r', label=r'$\sigma_3(t=0)$')
ax1.plot(s3.time_d, s3['S3 Magnitude']/1e6, c='k', alpha=0.5, lw=4, label=r"$\sigma_3$")
ax1.plot(sxx.time_d, sxx['Total Stress']/1e6, c='blue',lw=1.2, ls='--', label=r"$\sigma_t$")
ax1.plot(szz.time_d, szz['Total Stress']/1e6, c='green', lw=1.2, ls='-.', label=r"$\sigma_{zz}$")

## Caprock
dfxx = df_cap[df_cap.ten9_comp=="xx"]
s3 = df_cap.groupby('time_d')['S3 Magnitude'].max().reset_index()
sxx = dfxx.groupby('time_d')['Total Stress'].max().reset_index()
dfzz = df_cap[df_cap.ten9_comp=="zz"]
szz = dfzz.groupby('time_d')['Total Stress'].max().reset_index()
s3_0 = s3[s3.time_d<0]['S3 Magnitude'].values[0]

ax2.axhline(s3_0/1e6, ls=':', lw=2, alpha=0.3, c='r', label=r'$\sigma_3^0$')
ax2.plot(s3.time_d, s3['S3 Magnitude']/1e6, c='k', alpha=0.5, lw=4, label=r"$\sigma_3$")
ax2.plot(sxx.time_d, sxx['Total Stress']/1e6, c='blue', lw=1.2, ls='--', label=r"$\sigma_t$")
ax2.plot(szz.time_d, szz['Total Stress']/1e6, c='green', lw=1.2, ls='-.',label=r"$\sigma_{zz}$")

# Decorations
ax1.set_title("Stress 5m into the reservoir")
ax2.set_title("Stress 5m into the caprock")
for ax in [ax1,ax2] : 
    ax.set_xlabel("Time (days)")
    ax.set_xlim(0,365*2)

ax2.legend(loc='upper left', bbox_to_anchor=(1.02,1), fontsize=10, frameon=False)
ax1.set_ylabel("Stress (MPa)")
ax1.invert_yaxis()

fig.savefig("png/stress_by_time.png", dpi=500)
plt.show()
