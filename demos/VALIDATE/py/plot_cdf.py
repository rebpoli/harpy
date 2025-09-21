#!/usr/bin/env -S python -i
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.lines import Line2D
from matplotlib import gridspec
from scipy.interpolate import griddata

from netcdf import read_netcdf

# -------------------------------------------------------------------
# Files to process
# -------------------------------------------------------------------
files = [
    {"fname": "run/cdf/plane_xy.cd", "axes": (0, 1), "labels": ("X", "Y")},
    {"fname": "run/cdf/plane_yz.cd", "axes": (1, 2), "labels": ("Y", "Z")},
    {"fname": "run/cdf/plane_yz_well.cd", "axes": (1, 2), "labels": ("Y", "Z")},
]

# -------------------------------------------------------------------
# Helper: create regularly spaced grid for quiver
# -------------------------------------------------------------------
def regular_grid(X, Y, n_points=15):
    """Return indices of points on a regular grid over X, Y."""
    xi = np.linspace(X.min(), X.max(), n_points)
    yi = np.linspace(Y.min(), Y.max(), n_points)
    Xg, Yg = np.meshgrid(xi, yi)
    # Find nearest neighbor indices
    idx = [np.argmin((X - x)**2 + (Y - y)**2) for x, y in zip(Xg.ravel(), Yg.ravel())]
    return np.array(idx)

# -------------------------------------------------------------------
# Load datasets
# -------------------------------------------------------------------
datasets = []
nsteps = None
max_nsteps = 5

for f in files:
    ds = read_netcdf(f["fname"])
    ax0, ax1 = f["axes"]

    A = ds["Coord"][:, ax0].values
    B = ds["Coord"][:, ax1].values

    if nsteps is None:
        nsteps = min(max_nsteps, len(ds.time))
    time = ds.time.values[:nsteps]

    S1 = ds["S1"].values[:nsteps]
    S3 = ds["S3"].values[:nsteps]
    mag = ds["S3 Magnitude"].values[:nsteps]

    # Regularly spaced quiver
    idx = regular_grid(A, B, n_points=15)
    print(idx)

    datasets.append({
        "meta": f,
        "time": time,
        "A": A, "B": B,
        "A_vec": A[idx], "B_vec": B[idx],
        "S1_vec": S1[:, idx, :], "S3_vec": S3[:, idx, :],
        "mag": mag, "mag_vec": mag[:, idx],
        "grid": (np.linspace(A.min(), A.max(), 400),
                 np.linspace(B.min(), B.max(), 400)),
    })

# -------------------------------------------------------------------
# Global min/max for colorscale
# -------------------------------------------------------------------
vmin = min(d["mag"].min() for d in datasets)
vmax = max(d["mag"].max() for d in datasets)

# -------------------------------------------------------------------
# Figure setup (full space)
# -------------------------------------------------------------------
ncols = len(datasets)

from matplotlib import gridspec

fig = plt.figure(figsize=(6*ncols, 6), constrained_layout=True)
gs = gridspec.GridSpec(1, ncols + 1, figure=fig, width_ratios=[1]*ncols + [0.05])

axes = [fig.add_subplot(gs[0, i]) for i in range(ncols)]
cax = fig.add_subplot(gs[0, -1])

plots = []

for ax, d in zip(axes, datasets):
    A, B = d["A"], d["B"]
    Ag, Bg = np.meshgrid(*d["grid"])

    # Background
    Mag0 = griddata((A, B), d["mag"][0], (Ag, Bg), method="cubic")
    im = ax.imshow(
        Mag0, extent=[A.min(), A.max(), B.min(), B.max()],
        origin="lower", cmap="coolwarm", vmin=vmin, vmax=vmax,
        aspect="auto"
    )

    width_pixels = 2
    axes_width_pixels = ax.get_window_extent().width
    width_fraction = width_pixels / axes_width_pixels

    # Quivers
    ax0, ax1 = d["meta"]["axes"]
    quiv3 = ax.quiver(d["A_vec"], d["B_vec"], d["S3_vec"][0,:,ax0], d["S3_vec"][0,:,ax1],
                      angles="xy", scale_units="dots", scale=1/30,
                      color="green", width=width_fraction, pivot="middle", alpha=0.7,
                      headwidth=0, headlength=0, headaxislength=0)
    quiv1 = ax.quiver(d["A_vec"], d["B_vec"], d["S1_vec"][0,:,ax0], d["S1_vec"][0,:,ax1],
                      angles="xy", scale_units="dots", scale=1/30,
                      color="black", width=width_fraction, pivot="middle", alpha=0.7,
                      headwidth=0, headlength=0, headaxislength=0)

    ax.set_xlabel(d["meta"]["labels"][0])
    ax.set_ylabel(d["meta"]["labels"][1])
    title = ax.set_title(f"S1 & S3 on {d['meta']['labels'][0]}{d['meta']['labels'][1]} plane (time={d['time'][0]:,.0f})")

    plots.append({"im": im, "quiv1": quiv1, "quiv3": quiv3, "title": title, "data": d})

# Colorbar
plt.colorbar(plots[0]["im"], cax=cax, label="Temperature")

# Legend
legend_elements = [Line2D([0], [0], color="black", lw=2, label="S1"),
                   Line2D([0], [0], color="green", lw=2, label="S3")]
axes[0].legend(handles=legend_elements, loc="upper right")

# -------------------------------------------------------------------
# Animation update
# -------------------------------------------------------------------
def update(frame):
    updated = []
    for p in plots:
        d = p["data"]
        A, B = d["A"], d["B"]
        Ag, Bg = np.meshgrid(*d["grid"])

        Mag = griddata((A, B), d["mag"][frame], (Ag, Bg), method="linear")
        p["im"].set_data(Mag)

        ax0, ax1 = d["meta"]["axes"]
        p["quiv3"].set_UVC(d["S3_vec"][frame, :, ax0], d["S3_vec"][frame, :, ax1])
        p["quiv1"].set_UVC(d["S1_vec"][frame, :, ax0], d["S1_vec"][frame, :, ax1])

        p["title"].set_text(f"S1 & S3 on {d['meta']['labels'][0]}{d['meta']['labels'][1]} plane (time={d['time'][frame]:,.0f})")
        updated.extend([p["im"], p["quiv1"], p["quiv3"], p["title"]])
    return updated

anim = FuncAnimation(fig, update, frames=nsteps, interval=50, blit=False)
# anim.save("anim.mp4", dpi=200)

plt.show()


