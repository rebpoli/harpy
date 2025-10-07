#!/usr/bin/env -S python -i

import pandas as pd
import matplotlib.pyplot as plt
import os, time
import numpy as np
from matplotlib.animation import FuncAnimation

fn = "py/paper.mplstyle"
if os.path.isfile(fn) :
    plt.style.use('default')   ## reset!
    plt.style.use(fn)

from netcdf import read_netcdf
ds = read_netcdf("01-Transient/run/cdf/lin_0.cd")

# Find the point closest to the center of the top Y face
coords = ds['Coord'].values  # shape: (point_idx, vec3_comp)
# Extract x, y, z coordinates
x_coords = coords[:, 0]
y_coords = coords[:, 1]
z_coords = coords[:, 2]

# Find the top Y face (maximum Y value)
max_y = y_coords.max()
tolerance = (y_coords.max() - y_coords.min()) * 0.01  # 1% tolerance
top_face_mask = y_coords >= (max_y - tolerance)
top_face_indices = np.where(top_face_mask)[0]

print(f"Number of points on top Y face: {len(top_face_indices)}")
print(f"Top Y value: {max_y:.3f}")

# Find the center of the top face (in X-Z plane)
top_face_coords = coords[top_face_mask]
center_x = top_face_coords[:, 0].mean()
center_z = top_face_coords[:, 2].mean()

print(f"Center of top face: X={center_x:.3f}, Y={max_y:.3f}, Z={center_z:.3f}")

# Find the point on the top face closest to the center
distances_xz = np.sqrt((top_face_coords[:, 0] - center_x)**2 + 
                       (top_face_coords[:, 2] - center_z)**2)
closest_on_top_idx = top_face_indices[distances_xz.argmin()]

print(f"\nClosest point to center of top face:")
print(f"  Index: {closest_on_top_idx}")
print(f"  Coordinates: X={coords[closest_on_top_idx, 0]:.3f}, "
      f"Y={coords[closest_on_top_idx, 1]:.3f}, "
      f"Z={coords[closest_on_top_idx, 2]:.3f}")

center_idx = closest_on_top_idx

# Get displacement data at the center point
U = ds['U']  # shape: (time, point_idx, vec3_comp)
U_center = U[:, center_idx, :]  # shape: (time, vec3_comp)

# Get the index for Y component
y_idx = list(ds['vec3_comp'].values).index('y')

# Extract UY (Y-component of displacement)
UY = U_center[:, y_idx]

# Calculate Y_range (height of the model)
Y_range = y_coords.max() - y_coords.min()
print(f"\nModel Y range (height): {Y_range:.6f}")

# Calculate Y strain as UY / Y_range
Y_strain = UY / Y_range

# Get time values and convert from seconds to days
time_seconds = ds['time'].values
time_days = time_seconds / 86400  # Convert seconds to days (86400 s/day)

# Select one point every 12 hours (0.5 days)
interval_days = 0.1  # 12 hours in days
selected_indices = [0,1]  # Always include the first point

for i in range(1, len(time_days)):
    if time_days[i] - time_days[selected_indices[-1]] >= interval_days:
        selected_indices.append(i)

print(f"\nTotal time points: {len(time_seconds)}")
print(f"Selected time points (every 12h): {len(selected_indices)}")

# Overwrite with filtered data
time_days = time_days[selected_indices]
Y_strain = Y_strain[selected_indices]

# Create the plot
plt.figure(figsize=(12, 6))
plt.scatter(time_days, Y_strain, marker='x', s=50, color='#2E86AB')
plt.xlabel('Time (days)', fontsize=13)
plt.ylabel('Y Strain (UY / Y_range)', fontsize=13)
plt.title('Y-Strain at Top Face Center vs Time', fontsize=14, fontweight='bold')
plt.grid(True, alpha=0.3, linestyle='--')
plt.tight_layout()
plt.show()

# Optional: Print some statistics
print(f"\nTime range: {time_days[0]:.2f} to {time_days[-1]:.2f} days")
print(f"\nY Strain Statistics:")
print(f"  Min: {Y_strain.min():.6e}")
print(f"  Max: {Y_strain.max():.6e}")
print(f"  Mean: {Y_strain.mean():.6e}")
print(f"  Final value: {Y_strain[-1]:.6e}")
