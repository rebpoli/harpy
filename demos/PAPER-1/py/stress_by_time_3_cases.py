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

# # Load dataset
# filename = 'plane_yz.cd'
# filepath = f"run/cdf/{filename}"
# print(f"Loading {filepath}...")
# dataset = read_netcdf(filepath)


# z = dataset['Coord'].sel(vec3_comp= 'z').values

# # Caprock
# data_cap = dataset.isel(point_idx = z==-5)
# df_cap = data_cap[["Total Stress","S3 Magnitude"]].to_dataframe().reset_index()

# # Resrevoir
# data_res = dataset.isel(point_idx = z==5)
# df_res = data_res[["Total Stress","S3 Magnitude"]].to_dataframe().reset_index()

# # Save csv.
# df_res.to_csv("csv/df_res.csv")
# df_cap.to_csv("csv/df_cap.csv")

# Do some plotting
fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(16/2.54,6/2.64), sharey=True)

all_cases = [

        { 'dname' : "04-FULL-MODEL",         'legend' : 'Full model',   'c' : 'k', 'ls':'-' },
        { 'dname' : "05-NOCREEP",            'legend' : 'No creep',     'c' : 'r', 'ls':':'  },
        { 'dname' : "06-NOPRESSURESOLUTION", 'legend': 'No PS creep' , 'c' : 'b', 'ls':'-.' },
        { 'dname' : "08-ISOTHERMAL",         'legend': 'Isothermal' , 'c' : 'orange', 'ls':'--' }
]

for reg in all_cases :
    import pandas as pd
    dname = reg['dname']
    df_res = pd.read_csv(f"{dname}/csv/df_res.csv")
    df_cap = pd.read_csv(f"{dname}/csv/df_cap.csv")

    ## Time in days
    df_res["time_d"] = df_res.time/60/60/24
    df_cap["time_d"] = df_cap.time/60/60/24


    ## RES
    dfxx = df_res[(df_res.ten9_comp=="xx") | (df_res.ten9_comp=="yy")]
    s3 = df_res.groupby('time_d')['S3 Magnitude'].max().reset_index()
    sxx = dfxx.groupby('time_d')['Total Stress'].max().reset_index()
    dfzz = df_res[df_res.ten9_comp=="zz"]
    szz = dfzz.groupby('time_d')['Total Stress'].max().reset_index()
    s3_0 = s3[s3.time_d<0]['S3 Magnitude'].values[0]

    ax1.plot(s3.time_d, s3['S3 Magnitude']/1e6, c=reg['c'], ls=reg['ls'], lw=1, label=reg['legend'])

    ## Caprock
    dfxx = df_cap[df_cap.ten9_comp=="xx"]
    s3 = df_cap.groupby('time_d')['S3 Magnitude'].max().reset_index()
    sxx = dfxx.groupby('time_d')['Total Stress'].max().reset_index()
    dfzz = df_cap[df_cap.ten9_comp=="zz"]
    szz = dfzz.groupby('time_d')['Total Stress'].max().reset_index()
    s3_0 = s3[s3.time_d<0]['S3 Magnitude'].values[0]

    ax2.plot(s3.time_d, s3['S3 Magnitude']/1e6, c=reg['c'], ls=reg['ls'], lw=1, label=reg['legend'])

# Decorations
ax1.set_title(r"Minimum principal stress 5m into the reservoir")
ax2.set_title(r"Minimum principal stress 5m into the caprock")
for ax in [ax1,ax2] : 
    ax.set_xlabel("Time (days)")
    ax.set_xlim(0,365*2)

# leg=ax1.legend(loc='upper center', bbox_to_anchor=(0.5,0.98),
#            ncol=4, columnspacing=1, fontsize=7, frameon=True, edgecolor='k', fancybox=False, framealpha=1)
leg=ax2.legend(loc='upper left', bbox_to_anchor=(1.05,1),
           columnspacing=1, fontsize=7, frameon=True, edgecolor='k', fancybox=False, framealpha=1)
leg.get_frame().set_linewidth(0.5)
ax1.set_ylabel(r"Minimum principal stress, $\sigma_3$ (MPa)")
ax1.invert_yaxis()

fig.savefig("png/stress_by_time_with_and_without_creep.png", dpi=500)
plt.show()
