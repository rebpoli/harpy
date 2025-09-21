#!/usr/bin/env -S python -i

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter, FFMpegWriter
from timestr import format_time_duration

from netcdf import read_netcdf

#
#
def extract_vertical_profiles(y_coords, z_coords, s3_magnitude, num_lines=10):
    # Find Y coordinate bounds
    y_min, y_max = y_coords.min(), y_coords.max()

    # Create evenly spaced Y positions for vertical profiles
    y_line_positions = np.linspace(y_min, y_max, num_lines)

    profiles = []

    for i, y_pos in enumerate(y_line_positions):
        # Find points near this Y position
        y_tolerance = (y_max - y_min) / (num_lines * 1.5)  # Overlap for better sampling
        mask = np.abs(y_coords - y_pos) <= y_tolerance

        if np.any(mask):
            # Get points in this vertical slice
            z_slice = z_coords[mask]
            s3_mag_slice = s3_magnitude[mask]

            # Sort by Z coordinate (bottom to top)
            sort_idx = np.argsort(z_slice)
            z_sorted = z_slice[sort_idx]
            s3_mag_sorted = s3_mag_slice[sort_idx]

            profiles.append({
                'y_position': y_pos,
                'z_coords': z_sorted,
                's3_magnitude': s3_mag_sorted,
                'line_label': f'Y={y_pos:.1f}'
            })

    return profiles

#
#
def create_s3_vertical_line_animation(filename='plane_yz.cd', max_timesteps=None,
                                     num_lines=8, save_animation=True, show_animation=True):
    print(f"Creating S3 magnitude vs Z animation from {filename}")
    print(f"Number of vertical profiles: {num_lines}")

    # Load dataset
    filepath = f"run/cdf/{filename}"
    print(f"Loading {filepath}...")
    ds = read_netcdf(filepath)

    # Get time values
    time_values = ds.time.values
    print(f"Found {len(time_values)} time steps")

    # Determine timesteps to animate
    n_times = len(time_values)
    if max_timesteps is not None:
        n_times = min(max_timesteps, len(time_values))
        time_indices = np.linspace(0, len(time_values)-1, n_times, dtype=int)
        time_values_anim = time_values[time_indices]
    else:
        time_indices = np.arange(len(time_values))
        time_values_anim = time_values

    print(f"Animating {n_times} timesteps")

    # Get coordinates (YZ plane: Y=index 1, Z=index 2)
    coords = ds['Coord'].values
    y_coords = coords[:, 1]  # Y coordinates
    z_coords = coords[:, 2]  # Z coordinates

    y_min, y_max = y_coords.min(), y_coords.max()
    z_min, z_max = z_coords.min(), z_coords.max()

    print(f"Coordinate range: Y=[{y_min:.1f}, {y_max:.1f}], Z=[{z_min:.1f}, {z_max:.1f}]")

    # Process all timesteps to find global S3 magnitude range
    print("Analyzing S3 magnitude range...")
    all_s3_magnitudes = []

    for t_idx in [0, n_times//2, -1] if n_times > 2 else [0]:
        s3_magnitude = ds['S3 Magnitude'][time_indices[t_idx], :].values /1e6
        valid_mag = s3_magnitude[~np.isnan(s3_magnitude)]
        all_s3_magnitudes.extend(valid_mag)

    s3_min = np.percentile(all_s3_magnitudes, 1)
    s3_max = np.percentile(all_s3_magnitudes, 99)

    # Add 5 MPa padding to X-axis range
    padding_mpa = 5  
    s3_min_padded = s3_min - padding_mpa
    s3_max_padded = s3_max + padding_mpa

    print(f"S3 magnitude range: [{s3_min:.2e}, {s3_max:.2e}]")
    print(f"S3 magnitude range (with 5 MPa padding): [{s3_min_padded:.2e}, {s3_max_padded:.2e}]")

    # Process all timesteps
    print("Processing timesteps...")
    frame_data_list = []

    for frame_idx, time_idx in enumerate(time_indices):
        current_time = time_values_anim[frame_idx]

        # Get S3 magnitude at current time (directly from dataset)
        s3_magnitude = ds['S3 Magnitude'][time_idx, :].values /1e6

        # Extract vertical profiles
        profiles = extract_vertical_profiles(y_coords, z_coords, s3_magnitude, num_lines)

        frame_data_list.append({
            'time_idx': time_idx,
            'time_value': current_time,
            'profiles': profiles
        })

        if frame_idx % 5 == 0:
            print(f"  Processed frame {frame_idx+1}/{n_times}")

    # Create figure with more space for title
    fig, ax = plt.subplots(figsize=(6, 8))

    # Define colors for different vertical lines
    colors = plt.cm.tab10(np.linspace(0, 1, num_lines))

    # Animation update function
    def update_frame(frame_idx):
        frame_data = frame_data_list[frame_idx]
        current_time = frame_data['time_value']

        ax.clear()

        # Reference
        frame0 = frame_data_list[0];
        for i, profile in enumerate(frame0['profiles']):
            z_coords = profile['z_coords']
            s3_magnitude = profile['s3_magnitude']
            ax.plot(s3_magnitude, z_coords, color='k', alpha=0.5, linewidth=2, ls='--', label=profile['line_label'])

        # Plot each vertical profile
        for i, profile in enumerate(frame_data['profiles']):
            z_coords = profile['z_coords']
            s3_magnitude = profile['s3_magnitude']
            ax.plot(s3_magnitude, z_coords, color=colors[i], linewidth=2, label=profile['line_label'])


        # Styling
        ax.set_xlabel('S3 (MPa)', fontsize=12)
        ax.set_ylabel('Z (m)', fontsize=12)
        ax.set_title(f'Stress equilibrium after salt creep\n{format_time_duration(current_time)}', fontsize=14, pad=10)

        # Set consistent axis limits with 5 MPa padding on X-axis
        ax.set_xlim(s3_min_padded, s3_max_padded)
        ax.set_ylim(z_min, z_max)

        # Add grid
        ax.grid(True, alpha=0.3)

        # Add legend
#         ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)

        return [ax]

    # Create animation
    print("Creating animation...")
    anim = FuncAnimation(fig, update_frame, frames=n_times,
                        interval=500, blit=False, repeat=True)

    # Adjust layout to fit legend and leave space for title
    plt.subplots_adjust(top=0.9, right=0.95)

    # Save animation
    if save_animation:
        filename_out = f'S3_vertical_profiles_{n_times}steps_{num_lines}lines.gif'
        print(f"Saving {filename_out}...")
        anim.save(filename_out, writer=PillowWriter(fps=2))
        print("Animation saved!")

    # Show animation
    if show_animation:
        filename = f'S3_log.mp4'
        print(f"Saving {filename}...")

        # Use FFmpeg writer with high quality settings
        writer = FFMpegWriter(
            fps=1000//200,           # Frame rate based on interval
            metadata={},
            bitrate=8000,                 # High bitrate for quality
            extra_args=['-vcodec', 'libx264',   # H.264 codec
                       '-pix_fmt', 'yuv420p',   # Pixel format for compatibility
                       '-crf', '18',            # High quality (lower = better)
                       '-preset', 'slow']       # Better compression
        )
        anim.save(filename, writer=writer, dpi=150)  # High DPI for resolution
        print("Animation saved!")

    return fig, anim

#
#
#
if __name__ == "__main__":
    filename = 'plane_yz.cd'

    fig, anim = create_s3_vertical_line_animation(
        filename='plane_yz.cd',
        max_timesteps=10000,
        num_lines=1,
        save_animation=False,
        show_animation=True
    )
