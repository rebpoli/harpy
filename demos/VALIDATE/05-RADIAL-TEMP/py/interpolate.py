#!/usr/bin/env -S python

import os, time
import re
import numpy as np
import pandas as pd


df = pd.read_pickle( "db/temp.pkl" )
df = df.reset_index()

# Get sorted unique grid axes
t_vals = np.sort(df['t'].unique())
x_vals = np.sort(df['x'].unique())
z_vals = np.sort(df['z'].unique())

# Create 3D grid of temperature values
temp_grid = np.full((len(t_vals), len(x_vals), len(z_vals)), np.nan)

# Create lookup maps
t_idx = {v: i for i, v in enumerate(t_vals)}
x_idx = {v: i for i, v in enumerate(x_vals)}
z_idx = {v: i for i, v in enumerate(z_vals)}

for row in df.itertuples():
    ti = t_idx[row.t]
    xi = x_idx[row.x]
    zi = z_idx[row.z]
    temp_grid[ti, xi, zi] = row.TEMP

# Handle any missing values here if needed (e.g., fillna/interpolate)

from scipy.interpolate import RegularGridInterpolator
temp_at = RegularGridInterpolator(
    (t_vals, x_vals, z_vals), temp_grid,
    bounds_error=False, fill_value=np.nan
)




# # print (df)

# # points = df[['t', 'x', 'z']].values      # shape (N, 2)
# # values = df['TEMP'].values          # shape (N,)

# # # Build interpolator
# # print("Building interpolator ...")
# # from scipy.interpolate import LinearNDInterpolator
# # temp_at = LinearNDInterpolator(points, values)

# print("Interpolating ...")
# t = temp_at( (0, 0.053, 5205) )
# print(t)


# ## Try it
# import matplotlib.pyplot as plt
# import numpy as np
# from matplotlib.ticker import FuncFormatter

# # Choose resolution for output heatmap
# nx = 10000
# nz = 300

# # Fixed z range
# x_grid = np.linspace(x_vals.min(), 200, nx)
# z_grid = np.linspace(z_vals.min(), z_vals.max(), nz)
# X, Z = np.meshgrid(x_grid, z_grid)

# for t in t_vals:
#     print(f"t={t}")

#     t_grid = np.full_like(X, t)
#     query_points = np.stack((t_grid, X, Z), axis=-1).reshape(-1, 3)

#     # Interpolate
#     T_grid = temp_at(query_points).reshape(Z.shape)

#     # Plot
#     plt.figure(figsize=(6, 5))
#     plt.imshow(
#         T_grid,
#         extent=[x_grid.min(), x_grid.max(), z_grid.min(), z_grid.max()],
#         origin='lower',
#         aspect='auto',
#         cmap='coolwarm',
#         vmin=temp_grid.min(),  # consistent scale
#         vmax=temp_grid.max()
#     )
#     plt.gca().xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:.1f}"))
#     plt.gca().invert_yaxis()
#     plt.colorbar(label="Temperature")
#     plt.xlabel("R (m)")
#     plt.ylabel("TVDSS (m)")
#     plt.title(f"Temperature at t = {t:.1f}")
#     plt.tight_layout()
#     plt.savefig(f"png/heatmap_t{t:.1f}.png", dpi=150)
#     plt.close()



