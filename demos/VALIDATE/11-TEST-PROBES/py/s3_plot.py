#!/usr/bin/env -S python -i
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

ds = xr.open_dataset("run/cdf/plane_0.cd")

for v in [ 'Coord', 'S1', 'S2', 'S3' ] :
    lab = ds[v].attrs['components'].split() 
    ds = ds.assign_coords(vec3_comp=lab)


import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.interpolate import griddata

# Coordinates (static)
Y = ds["Coord"][:, 1].values
Z = ds["Coord"][:, 2].values

subsample = 4

# Subsample the 1D coordinates directly
Y_vec = Y[::subsample]
Z_vec = Z[::subsample]

# If you need a 2D meshgrid for quiver:
Y_grid, Z_grid = np.meshgrid(Y_vec, Z_vec, indexing='ij')

# --- Extract all timesteps ---
S1 = ds["S1"].values          # shape (time, point_idx, vec3_comp)
S3 = ds["S3"].values          # shape (time, point_idx, vec3_comp)
mag = ds["Delta_T"].values

# Subsample vectors
# Make sure point_idx corresponds to flattened grid order (row-major)
S1_vec = S1[:, ::subsample, :]
S3_vec = S3[:, ::subsample, :]

# --- Build grid for interpolation ---
ny, nz = 200, 200  # resolution of interpolated grid
y_lin = np.linspace(Y.min(), Y.max(), ny)
z_lin = np.linspace(Z.min(), Z.max(), nz)
Yg, Zg = np.meshgrid(y_lin, z_lin)

# --- Interpolate magnitude for first frame ---
Mag0 = griddata((Y, Z), mag[0], (Yg, Zg), method="linear")

# --- Setup figure ---
fig, ax = plt.subplots(figsize=(8, 6))

# Background interpolated colormap
im = ax.imshow(
    Mag0, extent=[Y.min(), Y.max(), Z.max(), Z.min()],
    origin="upper", cmap="coolwarm",
    vmin=mag.min(), vmax=mag.max(),
    aspect="auto"
)

# Vectors (first frame)
Vy = S3_vec[0, :, 1]
Vz = S3_vec[0, :, 2]
quiv = ax.quiver(
    Y_vec, Z_vec, Vy, Vz,
    angles="xy", scale_units="xy", scale=0.5,
    color="black", width=0.0015,
    pivot='middle', alpha=0.5,
    headwidth=0, headlength=0, headaxislength=0
)

Vy = S1_vec[0, :, 1]
Vz = S1_vec[0, :, 2]
quiv2 = ax.quiver(
    Y_vec, Z_vec, Vy, Vz,
    angles="xy", scale_units="xy", scale=0.5,
    color="white", width=0.0015,
    pivot='middle', alpha=0.5,
    headwidth=0, headlength=0, headaxislength=0
)

ax.set_xlabel("Y")
ax.set_ylabel("Z")
title = ax.set_title(f"S1&S3 projected on YZ plane (time={ds.time.values[0]})")
cbar = plt.colorbar(im, ax=ax, label="Temperatre")

# --- Update function ---
def update(frame):
    Mag = griddata((Y, Z), mag[frame], (Yg, Zg), method="linear")
    im.set_data(Mag)
    Vy = S3_vec[frame, :, 1]
    Vz = S3_vec[frame, :, 2]
    quiv.set_UVC(Vy, Vz)

    Vy = S1_vec[frame, :, 1]
    Vz = S1_vec[frame, :, 2]
    quiv2.set_UVC(Vy, Vz)
    title.set_text(f"S1&S3 projected on YZ plane (time={ds.time.values[frame]:,.0f})")
    return im, quiv, title

# --- Animate ---
anim = FuncAnimation(fig, update, frames=len(ds.time), interval=500, blit=False)

plt.show()

