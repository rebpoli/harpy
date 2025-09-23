#!/usr/bin/env python3

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
    find_closest_xy_positions, extract_stress_data_at_time, extract_vertical_profile,
    extract_time_series_at_position, find_closest_profile_for_depth, calculate_stress_ranges,
    create_profile_legend_elements, calculate_stress_differences, create_video_with_ffmpeg
)

# Set matplotlib backend for parallel processing
import matplotlib
matplotlib.use('Agg')

# Font settings
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.default'] = 'regular'
plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['figure.titlesize'] = 20


# ==================== GLOBAL CONTROL VARIABLES ====================
NUM_LINES = 5  # Number of profile lines to display
TRACKED_DEPTHS = [5.0, -5.0]  # Depths to track (reservoir, caprock)
WELL_POSITION = (0.0, 0.0)  # Well position (x, y)
DEPTH_COLORS = ['red', 'blue']  # Colors for tracked depths
STRESS_COLORS = ['purple', 'orange']  # Colors for sigmaxx, sigmazz
XY_TOLERANCE = 1e-6  # Tolerance for (x,y) position matching
Z_TOLERANCE = 2.0  # Tolerance for depth matching (meters)

#
#
#
def plot_vertical_profiles(ax, dataset, time_idx, xy_positions, time_series_data, stress_ranges):
    coords = dataset['Coord'].values
    z_coords = coords[:, 2]
    z_min, z_max = z_coords.min(), z_coords.max()
    sxx_min, sxx_max = stress_ranges

    # Extract current and reference stress data
    current_stress = extract_stress_data_at_time(dataset, time_idx)
    reference_stress = extract_stress_data_at_time(dataset, 0)  # First timestep

    # Plot reference profiles (gray)
    reference_profiles = []
    for xy_pos in xy_positions:
        profile = extract_vertical_profile(reference_stress, xy_pos)
        if profile is not None:
            reference_profiles.append(profile)
            z_coords_prof = profile['z_coords']
            sxx_prof = profile['sxx_magnitude']
            sigmazz_prof = profile['sigmazz_magnitude']

            ax.plot(sxx_prof, z_coords_prof, color='orange', alpha=0.5, linewidth=1, ls='-')
            ax.plot(sigmazz_prof, z_coords_prof, color='orange', alpha=0.5, linewidth=1, ls='--')

    # Plot current profiles (colored)
    current_profiles = []
    for xy_pos in xy_positions:
        profile = extract_vertical_profile(current_stress, xy_pos)
        if profile is not None:
            current_profiles.append(profile)
            z_coords_prof = profile['z_coords']
            sxx_prof = profile['sxx_magnitude']
            sigmazz_prof = profile['sigmazz_magnitude']

            # Get color and intensity from xy_position
            color = xy_pos['color']
            intensity = xy_pos['intensity']
            alpha = 0.1 + 0.9 * intensity
            linewidth = 1

            # Plot sigmaxx (solid) and sigmazz (dashed)
            ax.plot(sxx_prof, z_coords_prof, color=color, linewidth=linewidth, alpha=alpha)
            ax.plot(sigmazz_prof, z_coords_prof, color=color, linewidth=linewidth, alpha=alpha, linestyle='--')

    # Add tracking markers
    for j, (target_z, depth_color) in enumerate(zip(TRACKED_DEPTHS, DEPTH_COLORS)):
        if target_z in time_series_data:
            # Use time series data for marker values
            sxx_current = time_series_data[target_z]['sigmaxx'][time_idx]
            sigmazz_current = time_series_data[target_z]['sigmazz'][time_idx]

            # sigmaxx circles
            ax.plot(sxx_current, target_z, 'o', markersize=12,
                   color=depth_color, markerfacecolor='white',
                   markeredgewidth=3, markeredgecolor=depth_color)

            # Sigmazz squares
            ax.plot(sigmazz_current, target_z, 's', markersize=10,
                   color=depth_color, markerfacecolor='white',
                   markeredgewidth=3, markeredgecolor=depth_color)

    # Configure axis
    ax.set_ylabel('Z (m)', fontsize=12)
    ax.set_xlim(sxx_max, sxx_min)
    ax.set_ylim(z_max, z_min)  # Inverted Z axis
    ax.grid(True, alpha=0.3)
    ax.tick_params(axis='x', labelbottom=False)

    # Add legend
    legend_elements = create_profile_legend_elements(TRACKED_DEPTHS)
    ax.legend(handles=legend_elements, fontsize=12, loc='upper left', ncol=1)

    return current_profiles

#
#
#
def plot_time_series(ax, time_series_data, time_years, current_time_years, time_idx, stress_ranges):
    sxx_min, sxx_max = stress_ranges

    # Plot time series for each depth
    for j, (target_z, depth_color) in enumerate(zip(TRACKED_DEPTHS, DEPTH_COLORS)):
        if target_z in time_series_data:
            depth_data = time_series_data[target_z]
            if target_z > 0 :
                depth_label = 'Reservoir'
                color = 'red'
            else:
                depth_label = 'Caprock'
                color = 'blue'

            ax.plot(depth_data['sigmaxx'], time_years,
                   color=color, linestyle='-', linewidth=1, alpha=0.8,
                   label =  r'$\sigma_{xx}$ ' + f'{depth_label}')

            ax.plot(depth_data['sigmazz'], time_years,
                   color=color, linestyle='--', linewidth=1, alpha=0.8,
                   label =  r'$\sigma_{zz}$ ' + f'{depth_label}')

            # Add current time markers
            current_sxx = depth_data['sigmaxx'][time_idx]
            current_sigmazz = depth_data['sigmazz'][time_idx]

            ax.plot(current_sxx, current_time_years, 'o',
                   markersize=10, color=depth_color,
                   markerfacecolor='white', markeredgewidth=3,
                   markeredgecolor=depth_color)

            ax.plot(current_sigmazz, current_time_years, 's',
                   markersize=8, color=depth_color,
                   markerfacecolor='white', markeredgewidth=3,
                   markeredgecolor=depth_color)

    # Configure axis
    ax.set_ylabel('Time (years)', fontsize=12)
    ax.set_xlabel(r'Stress (MPa)', fontsize=12)
#     print(f"YLIM: {time_years[-1]} - {time_years[0]}")
#     ax.set_ylim(1, 0)  # Invert time axis - 1 year only
    ax.set_ylim(time_years[-1], time_years[0])  # Invert time axis
    ax.set_xlim(sxx_max, sxx_min)
    ax.grid(True, alpha=0.3)

    # Add horizontal line at current time
    ax.axhline(y=current_time_years, color='gray', linestyle='-', alpha=0.8, linewidth=1)

    # Add legend
    ax.legend(loc='lower left', fontsize=12)

    # Add stress difference text
    stress_diffs = calculate_stress_differences(time_series_data, TRACKED_DEPTHS, time_idx)
    if stress_diffs is not None:
        line1 = r'$\Delta\sigma$ Cap-Res (hh): ' + f'{stress_diffs["delta_between"]:+.1f} MPa'
        line2 = r'$\Delta\sigma$ Caprock (hv): ' + f'{stress_diffs["delta_caprock"]:+.1f} MPa'
        line3 = r'$\Delta\sigma$ Reservoir (hv): ' + f'{stress_diffs["delta_reservoir"]:+.1f} MPa'

        diff_text = f'{line1}\n{line2}\n{line3}'

        text_x = (sxx_min + sxx_max) / 2
        dtmax = time_years[-1] - time_years[0]
        text_y = time_years[-1] - dtmax*0.05 #time_years[0] + dtmax * 0.4

        ax.text(text_x, text_y, diff_text,
               fontsize=12, ha='center', va='bottom',
               bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9))

#
#
#
def process_timestep_and_save_frame(args):
    (frame_idx, time_idx, dataset, time_values_anim, time_years_anim,
     xy_positions, time_series_data, stress_ranges, output_dir, dpi) = args
    coords = dataset['Coord'].values

    # Extract current time information
    current_time = time_values_anim[frame_idx]
    current_time_years = time_years_anim[frame_idx]
    time_years = time_years_anim

    # Create figure
    width_px, height_px = 1600, 2050
    dpi = 150
    figsize = (width_px/dpi, height_px/dpi)  # (19.2, 10.8) inches
    fig, (ax_profiles, ax_timeseries) = plt.subplots(2, 1, figsize=figsize, dpi=dpi)

    # Add figure title
    fig.suptitle(f'Time: {format_time_duration(current_time)}',
                fontsize=16, x=0.02, y=0.98, ha='left', va='top')

    # Plot vertical profiles
    plot_vertical_profiles(ax_profiles, dataset, time_idx, xy_positions,
                          time_series_data, stress_ranges)

    # Plot time series
    plot_time_series(ax_timeseries, time_series_data, time_years,
                    current_time_years, time_idx, stress_ranges)

    # Adjust layout and save
    plt.subplots_adjust(left=0.12, right=0.95, bottom=0.08, top=0.95, hspace=0.05)

    frame_filename = os.path.join(output_dir, f'frame_{frame_idx:04d}.png')
    fig.savefig(frame_filename, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close(fig)

    return frame_filename

#
#  The main function
#
def create_combined_animation_parallel(dataset, max_timesteps=None, interval=300, dpi=150):
    print(f"Creating combined Sigma XX profiles and time series animation from {filename}")
    print(f"Number of vertical profiles: {NUM_LINES}")

    # Extract basic information
    time_values = dataset.time.values
    time_years = time_values / (365.25 * 24 * 3600)

    print(f"Found {len(time_values)} time steps")

    # Determine timesteps to animate
    n_times = len(time_values)
    if max_timesteps is not None:
        n_times = min(max_timesteps, len(time_values))
        time_indices = np.linspace(0, len(time_values)-1, n_times, dtype=int)
        time_values_anim = time_values[time_indices]
        time_years_anim = time_years[time_indices]
    else:
        time_indices = np.arange(len(time_values))
        time_values_anim = time_values
        time_years_anim = time_years

    print(f"Animating {n_times} timesteps")

    # Find closest (x,y) positions to well
    print("Finding closest positions to well...")
    xy_positions = find_closest_xy_positions(dataset, WELL_POSITION, NUM_LINES)

    for i, pos in enumerate(xy_positions):
        print(f"  Position {i}: ({pos['x_pos']:.2f}, {pos['y_pos']:.2f}), distance: {pos['distance']:.2f}m")

    # Extract time series data for the closest position
    print("Computing time series data for closest position to well...")
    closest_position = xy_positions[0]  # Use closest position for time series
    time_series_data = extract_time_series_at_position(dataset, closest_position, TRACKED_DEPTHS)
    print(f"Closest position: {closest_position}")

    print(f"Time series data extracted for {len(time_series_data)} depths")
    for depth in time_series_data:
        data = time_series_data[depth]
        print(f"  Depth {depth}m: Sigma XX [{np.min(data['sigmaxx']):.1f}, {np.max(data['sigmaxx']):.1f}] MPa")

    # Calculate stress ranges
    stress_ranges = calculate_stress_ranges(time_series_data)
    print(f"Stress range: [{stress_ranges[0]:.1f}, {stress_ranges[1]:.1f}] MPa")

    # Create output directory
    output_dir = 'temp_frames_sxx'
    os.makedirs(output_dir, exist_ok=True)

    # Prepare arguments for parallel processing
    print(f"Processing timesteps and generating frames in parallel...")
    process_args = [
        (frame_idx, time_indices[frame_idx], dataset, time_values_anim, time_years_anim,
         xy_positions, time_series_data, stress_ranges, output_dir, dpi)
        for frame_idx in range(n_times)
    ]

    # Process frames in parallel
    max_workers = min(cpu_count() - 1, 12)
    completed = 0
    frame_files = [None] * n_times
    start_time = time.time()

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_to_idx = {
            executor.submit(process_timestep_and_save_frame, args): args[0]
            for args in process_args
        }

        for future in as_completed(future_to_idx):
            frame_idx = future_to_idx[future]
            try:
                frame_file = future.result()
                frame_files[frame_idx] = frame_file
                completed += 1

                progress_percent = (completed / n_times) * 100
                elapsed_time = time.time() - start_time

                print(f"\rProcessing: {completed}/{n_times} ({progress_percent:.1f}%)", end="", flush=True)

            except Exception as exc:
                print(f"\nFrame {frame_idx} generated an exception: {exc}")

    print(f"\nGenerated {completed} frames in {time.time() - start_time:.1f}s")

    # Create video using FFmpeg
    fps = 1000 // interval
    output_filename = 'combined_stress_animation.mp4'
    success = create_video_with_ffmpeg(output_dir, output_filename, fps, cleanup_frames=True)
    
    if success:
        print(f"Animation saved as {output_filename}")
        return output_filename
    else:
        print("Animation creation failed.")
        return None

#
#
#
if __name__ == "__main__":
    # Load dataset
    filename = 'plane_yz.cd'
    filepath = f"run/cdf/{filename}"
    print(f"Loading {filepath}...")
    dataset = read_netcdf(filepath)

    # Work
    video_file = create_combined_animation_parallel(
        dataset,
        max_timesteps=1000, interval=200, dpi=150
    )

    print("Done")
