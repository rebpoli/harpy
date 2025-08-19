#!/usr/bin/env -S python3 -i

import argparse
import shutil
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator, ScalarFormatter
from matplotlib import colors as mcolors
from concurrent.futures import ProcessPoolExecutor
import subprocess

# ---------------------------------------------------------------------
#                         CONFIGURABLE PARAMETERS
# ---------------------------------------------------------------------
phi_salt_deg = 40.0       # Friction angle for salt [degrees]
phi_carb_deg = 50.0       # Friction angle for carbonate [degrees]
cohesion_salt = 2.0       # Cohesion for salt [MPa]
cohesion_carb = 10.0      # Cohesion for carbonate [MPa]
p_envelope_max = 80.0     # Maximum p' to draw the envelope up to [MPa]

# Convert to radians once
phi_salt = np.radians(phi_salt_deg)
phi_carb = np.radians(phi_carb_deg)

# Paper-like style
plt.style.use("seaborn-v0_8-paper")


# ---------------------------------------------------------------------
#                         HELPER FUNCTIONS
# ---------------------------------------------------------------------
def format_time(days: float) -> str:
    """Convert time in days to 'Xy Yd HH:MM'."""
    total_minutes = int(round(days * 24 * 60))
    years, rem = divmod(total_minutes, int(365.25 * 24 * 60))
    days_i, rem = divmod(rem, 24 * 60)
    hours, minutes = divmod(rem, 60)
    parts = []
    if years > 0:
        parts.append(f"{years}y")
    if days_i > 0:
        parts.append(f"{days_i}d")
    parts.append(f"{hours:02}:{minutes:02}")
    return " ".join(parts)


def style_axes(ax):
    ax.grid(True, which="major", linestyle="--", linewidth=0.6, alpha=0.85)
    ax.grid(True, which="minor", linestyle=":", linewidth=0.5, alpha=0.55)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.tick_params(axis="both", which="major", labelsize=11, width=1.2, length=6)
    ax.tick_params(axis="both", which="minor", width=1.0, length=3)
    for s in ax.spines.values():
        s.set_linewidth(1.2)


def mohr_coulomb_envelope(phi, cohesion, pmax, npoints=200):
    """
    Compute Mohr–Coulomb envelope in p-q space:
        q = p * tan(phi) + c / cos(phi)
    """
    p = np.linspace(0, pmax, npoints)
    q = p * np.tan(phi) + cohesion / np.cos(phi)
    return p, q


def build_dataframe_from_csv(csv_path: Path) -> pd.DataFrame:
    """Read the TSV and return a ready-to-plot DataFrame with columns: p,q,t,dz,z."""
    print(f"[1/7] Reading CSV: {csv_path}")
    df_src = pd.read_csv(csv_path, sep="\t", engine="python", quotechar='"')

    pressure = pd.to_numeric(df_src["Pressure"], errors="coerce").to_numpy()
    invar_p  = pd.to_numeric(df_src["Invar P"], errors="coerce").to_numpy()
    q_raw    = pd.to_numeric(df_src["Invar Q"], errors="coerce").to_numpy()
    t        = pd.to_numeric(df_src["Time(day)"], errors="coerce").to_numpy()
    z        = pd.to_numeric(df_src["Z"], errors="coerce").to_numpy()

    p = (pressure + invar_p) / 1e6  # MPa
    q = q_raw / 1e6                 # MPa
    dz = np.abs(z)

    df = pd.DataFrame({
        "p": p,
        "q": q,
        "t": t,
        "dz": dz,
        "z": z,
    })
    print(f"[2/7] Built DataFrame with {len(df)} rows.")
    return df


def load_or_build_df(csv_path: Path, pkl_path: Path, force: bool) -> pd.DataFrame:
    """Load cached pickle if fresh; otherwise rebuild from CSV and save."""
    csv_mtime = csv_path.stat().st_mtime
    need_rebuild = True

    if pkl_path.exists() and not force:
        pkl_mtime = pkl_path.stat().st_mtime
        if pkl_mtime >= csv_mtime:
            print(f"[Cache] Using up-to-date pickle: {pkl_path}")
            try:
                return pd.read_pickle(pkl_path)
            except Exception as e:
                print(f"[Cache] Failed to read pickle ({e}); rebuilding from CSV...")
        else:
            print(f"[Cache] Pickle older than CSV; rebuilding...")
    elif force:
        print(f"[Cache] --force-pkl requested; rebuilding from CSV...")
    else:
        print(f"[Cache] No pickle found; building from CSV...")

    df = build_dataframe_from_csv(csv_path)
    try:
        df.to_pickle(pkl_path)
        print(f"[Cache] Saved fresh pickle: {pkl_path}")
    except Exception as e:
        print(f"[Cache] Warning: failed to save pickle ({e}). Continuing without cache.")
    return df


def plot_frame(args):
    (
        ts, pvals, qvals, dzvals, zvals,
        xmin, xmax, ymin, ymax,
        dzmax, dpi, outdir
    ) = args

    carbonate = zvals > 0
    salt = ~carbonate

    norm = mcolors.Normalize(vmin=0.0, vmax=dzmax)
    cmap = plt.colormaps["coolwarm_r"]
    edge_salt = cmap(norm(np.clip(dzvals[salt], 0.0, dzmax)))
    edge_carb = cmap(norm(np.clip(dzvals[carbonate], 0.0, dzmax)))

    # Figure layout: 1 row × 3 cols (last for colorbar)
    fig = plt.figure(figsize=(10.5, 5.5), constrained_layout=False)
    gs = gridspec.GridSpec(nrows=1, ncols=3, width_ratios=[1, 1, 0.05], figure=fig)

    ax1 = fig.add_subplot(gs[0, 0])      # Salt
    ax2 = fig.add_subplot(gs[0, 1])      # Carbonate
    cax = fig.add_subplot(gs[0, 2])      # Colorbar

    # Salt points
    ax1.scatter(
        pvals[salt], qvals[salt],
        s=14, facecolors="none", edgecolors=edge_salt,
        marker="o", linewidths=0.9, alpha=0.95
    )
    ax1.set_title("Salt (Z ≤ 0)", fontsize=13)

    # Carbonate points
    ax2.scatter(
        pvals[carbonate], qvals[carbonate],
        s=14, facecolors="none", edgecolors=edge_carb,
        marker="o", linewidths=0.9, alpha=0.95
    )
    ax2.set_title("Carbonate (Z > 0)", fontsize=13)

    # Add Mohr–Coulomb envelopes
    p_env_salt, q_env_salt = mohr_coulomb_envelope(phi_salt, cohesion_salt, p_envelope_max)
    p_env_carb, q_env_carb = mohr_coulomb_envelope(phi_carb, cohesion_carb, p_envelope_max)

    for ax, (p_env, q_env), color in zip(
        (ax1, ax2),
        ((p_env_salt, q_env_salt), (p_env_carb, q_env_carb)),
        ("#1f77b4", "#ff7f0e")
    ):
        ax.plot(p_env, q_env, color=color, linestyle="--", linewidth=1.5, label="MC envelope")
        ax.legend(fontsize=9, loc="upper left", frameon=True)

    # Shared axes
    for ax in (ax1, ax2):
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax) # CHANGED ymax or xmax for a square scale
        ax.set_xlabel("p' (MPa)", fontsize=13)
        style_axes(ax)
    ax1.set_ylabel("q (MPa)", fontsize=13)

    # Title
    fig.suptitle(f"p–q diagram\nTime = {format_time(ts)}", fontsize=14)

    # Colorbar
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax, orientation="vertical")
    cbar.set_label("Distance to salt–carbonate interface (m)", fontsize=11)
    cbar.ax.tick_params(labelsize=10)

    # Spacing
    fig.subplots_adjust(left=0.08, right=0.97, top=0.86, bottom=0.12, wspace=0.22)
    pos = cbar.ax.get_position()
    cbar.ax.set_position([
        pos.x0 - 0.03, pos.y0,
        pos.width, pos.height
    ])

    fname = outdir / f"frame_{ts:010.5f}.png"
    plt.savefig(fname, dpi=dpi)
    plt.close(fig)
    return fname


# ---------------------------------------------------------------------
#                         CLI ENTRY POINT
# ---------------------------------------------------------------------
ap = argparse.ArgumentParser(
    description="Parallel p–q animation (Salt vs Carbonate; hollow markers; ParaView colors; cached DataFrame + MC envelopes)."
)
ap.add_argument("csv", type=Path, help="Path to TSV")
ap.add_argument("--output", type=Path, default=None, help="Output video path")
ap.add_argument("--fps", type=int, default=10, help="Frames per second")
ap.add_argument("--dpi", type=int, default=120, help="PNG DPI")
ap.add_argument("--crf", type=int, default=18, help="Quality (NVENC cq / VP9 crf)")
ap.add_argument("--skip", type=int, default=1, help="Use every Nth unique time to speed up")
ap.add_argument("--workers", type=int, default=20, help="Parallel workers for frame rendering")
ap.add_argument("--pkl", type=Path, default=None, help="Path to cache pickle (defaults next to CSV)")
ap.add_argument("--force-pkl", dest="force_pkl", action="store_true", help="Rebuild the pickle from CSV")
ap.add_argument("--force_pkl", dest="force_pkl", action="store_true", help="Rebuild the pickle from CSV (alias)")
args = ap.parse_args()

# Cache / load dataframe
pkl_path = args.pkl or args.csv.with_suffix("").with_name(args.csv.stem + "_prepared.pkl")
df = load_or_build_df(args.csv, pkl_path, force=args.force_pkl)
print(f"[3/7] DataFrame shape: {df.shape}  (cache: {pkl_path})")

# Group by time
groups = df.groupby("t")
unique_times = np.array(sorted(groups.groups.keys()), dtype=float)
if args.skip > 1:
    unique_times = unique_times[::args.skip]
nframes = len(unique_times)
print(f"[4/7] Using {nframes} frames (skip={args.skip}).")

xmin, xmax = 0.0, (float(df["p"].max()) * 1.05) if not df.empty else 1.0
ymin, ymax = 0.0, (float(df["q"].max()) * 1.05) if not df.empty else 1.0
dzmax = float(np.percentile(df["dz"], 99.0)) if len(df) > 0 else 1.0

outdir = Path("_frames")
if outdir.exists():
    shutil.rmtree(outdir)
outdir.mkdir(parents=True)

frame_args = []
for ts in unique_times:
    g = groups.get_group(ts)
    frame_args.append((
        ts,
        g["p"].to_numpy(),
        g["q"].to_numpy(),
        g["dz"].to_numpy(),
        g["z"].to_numpy(),
        xmin, xmax, ymin, ymax,
        dzmax,
        args.dpi, outdir
    ))

# Render frames
print(f"[5/7] Rendering {nframes} frames in parallel using {args.workers} workers...")
with ProcessPoolExecutor(max_workers=args.workers) as ex:
    list(ex.map(plot_frame, frame_args))

# Encode video
out_path = args.output or args.csv.with_suffix("").with_name(args.csv.stem + "_pq_anim.mp4")
ffmpeg_cmd = [
    "ffmpeg", "-y",
    "-framerate", str(args.fps),
    "-pattern_type", "glob",
    "-i", str(outdir / "*.png"),
    "-c:v", "h264_nvenc",
    "-pix_fmt", "yuv420p",
    "-cq", str(args.crf),
    "-b:v", "0",
    "-gpu", "0",
    str(out_path)
]

print(f"[6/7] Encoding video using FFmpeg (h264_nvenc) ...")
subprocess.run(ffmpeg_cmd, check=True)

print(f"[7/7] ✅ Saved animation to: {out_path}")
shutil.rmtree(outdir)
print("🧹 Temporary frames cleaned up.")
