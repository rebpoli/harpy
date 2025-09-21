#!/usr/bin/env -S python -i
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.interpolate import griddata

from subsample import spatial_subsample
from netcdf import read_netcdf
ds = read_netcdf( "run/cdf/plane_xy.cd" )

# Coordinates (static)
Y = ds["Coord"][:, 1].values
Z = ds["Coord"][:, 0].values


min_dist = 0.04 * max(Y.max()-Y.min(), Z.max()-Z.min())  # ~5% of domain size
idx = spatial_subsample(Y, Z, min_dist)

Y_vec = Y[idx]
Z_vec = Z[idx]

# If needed for quiver grid:
Y_grid, Z_grid = np.meshgrid(Y_vec, Z_vec, indexing='ij')

nsteps = min(80, len(ds.time))
time = ds.time.values[:nsteps]
S1 = ds["S1"].values[:nsteps]
S3 = ds["S3"].values[:nsteps]
# mag = ds["Delta_T"].values[:nsteps]
# mag = ds["Invariant P(eff)"].values[:nsteps]
mag = ds["S3 Magnitude"].values[:nsteps]

# Apply same selection to fields
S1_vec = S1[:, idx, :]
S3_vec = S3[:, idx, :]
mag_vec = mag[:, idx]


# Background grid
ny, nz = 400, 400
y_lin = np.linspace(Y.min(), Y.max(), ny)
z_lin = np.linspace(Z.min(), Z.max(), nz)
Yg, Zg = np.meshgrid(y_lin, z_lin)

# First frame background
Mag0 = griddata((Y, Z), mag[0], (Yg, Zg), method="cubic")

# --- Setup figure ---
fig, ax = plt.subplots(figsize=(8, 6))

# Background
im = ax.imshow(
    Mag0, extent=[Y.min(), Y.max(), Z.max(), Z.min()],
    origin="upper", cmap="coolwarm",
    vmin=mag.min(), vmax=mag.max(),
    aspect="auto"
)

# --- Quiver setup (S3: white, S1: black semi-transparent) ---
Vy3, Vz3 = S3_vec[0, :, 1], S3_vec[0, :, 0]
quiv3 = ax.quiver(
    Y_vec, Z_vec, Vy3, Vz3,
    angles="xy", scale_units="xy", scale=0.6,  # auto-scaling
    color="white", width=0.002,
    pivot="middle", alpha=0.7,
    headwidth=0, headlength=0, headaxislength=0  # tiny arrows
)

Vy1, Vz1 = S1_vec[0, :, 1], S1_vec[0, :, 0]
quiv1 = ax.quiver(
    Y_vec, Z_vec, Vy1, Vz1,
    angles="xy", scale_units="xy", scale=0.3,
    color="black", width=0.004,
    pivot="middle", alpha=0.7,
    headwidth=0, headlength=0, headaxislength=0
)

from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], color="black", lw=2, label="S1"),
    Line2D([0], [0], color="white", lw=2, label="S3")
]

ax.legend(handles=legend_elements, loc="upper right")

ax.set_xlabel("Y")
ax.set_ylabel("Z")
title = ax.set_title(f"S1 & S3 projected on YZ plane (time={ds.time.values[0]:,.0f})")
cbar = plt.colorbar(im, ax=ax, label="Temperature")

# --- Update function ---
def update(frame):
    Mag = griddata((Y, Z), mag[frame], (Yg, Zg), method="linear")
    im.set_data(Mag)

    quiv3.set_UVC(S3_vec[frame, :, 1], S3_vec[frame, :, 0])
    quiv1.set_UVC(S1_vec[frame, :, 1], S1_vec[frame, :, 0])

    title.set_text(f"S1 & S3 projected on YZ plane (time={ds.time.values[frame]:,.0f})")
    return im, quiv3, quiv1, title

# Animate
anim = FuncAnimation(fig, update, frames=nsteps, interval=10, blit=False)

# anim.save( "s1_s3_animation.mp4", writer="ffmpeg", fps=2, dpi=150, bitrate=-1)

plt.show()

