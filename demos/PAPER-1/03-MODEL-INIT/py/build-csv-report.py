#!/usr/bin/env -S python -i

# import os
# import subprocess
# import shutil
# import numpy as np
# import matplotlib
# import matplotlib.pyplot as plt
# from matplotlib.lines import Line2D
# from scipy.interpolate import griddata
# from concurrent.futures import ProcessPoolExecutor, as_completed
# from multiprocessing import cpu_count
# import time
# from plot_util import savefig

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from netcdf import read_netcdf
ds = read_netcdf("run/cdf/plane_yz.cd")

x_coord = ds['Coord'].isel(vec3_comp=0).values
y_coord = ds['Coord'].isel(vec3_comp=1).values
z_coord = ds['Coord'].isel(vec3_comp=2).values

# Find the middle X and Y
x_middle = (x_coord.max() + x_coord.min()) / 2
y_middle = (y_coord.max() + y_coord.min()) / 2
distance_from_middle = np.sqrt((x_coord - x_middle)**2 + (y_coord - y_middle)**2)

# Find unique Z values and select closest (X,Y) for each Z
unique_z = np.unique(z_coord)
selected_indices = []

for z_val in unique_z:
    # Find all points with this Z value
    z_mask = z_coord == z_val
    z_indices = np.where(z_mask)[0]

    # Find the point closest to middle among these
    closest_idx = z_indices[np.argmin(distance_from_middle[z_indices])]
    selected_indices.append(closest_idx)

selected_indices = np.array(selected_indices)

# Extract Total Stress components for selected points
total_stress = ds['Total Stress'].isel(time=-1, point_idx=selected_indices)
stress_xx = total_stress.isel(ten9_comp=0).values
stress_yy = total_stress.isel(ten9_comp=4).values
stress_zz = total_stress.isel(ten9_comp=8).values

pressure = ds['Pressure'].isel(time=-1, point_idx=selected_indices).values
temperature = ds['Temperature'].isel(time=-1, point_idx=selected_indices).values
delta_p = ds['Delta_P'].isel(time=-1, point_idx=selected_indices).values
delta_t = ds['Delta_T'].isel(time=-1, point_idx=selected_indices).values

# Create the dataframe
df = pd.DataFrame({
    'Z': z_coord[selected_indices],
    'X': x_coord[selected_indices],
    'Y': y_coord[selected_indices],
    'sigma_xx': stress_xx,
    'sigma_yy': stress_yy,
    'sigma_zz': stress_zz,
    'Pressure': pressure,
    'Temperature': temperature,
    'Delta_P': delta_p,
    'Delta_T': delta_t
})

df.to_pickle("pkl/stress_profile.pkl")
