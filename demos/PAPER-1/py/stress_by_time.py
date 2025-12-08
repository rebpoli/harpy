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

# Set matplotlib backend for parallel processing
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

fig, ax = plt.subplots(1, 1, figsize=(8/2.54,8/2.64))

z = dataset['Coord'].sel(vec3_comp= 'z').values

# Caprock
data_cap = dataset.isel(point_idx = z==-5)
df_cap = data_cap[["Total Stress","S3 Magnitude"]].to_dataframe().reset_index()

# Resrevoir
data_res = dataset.isel(point_idx = z==5)
df_res = data_cap[["Total Stress","S3 Magnitude"]].to_dataframe().reset_index()

# Save csv.
df_res.to_csv("csv/df_res.csv")
df_cap.to_csv("csv/df_cap.csv")

# Do some plotting
