#!/usr/bin/env python3
"""
Helper functions for Stress stress analysis and animation
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import subprocess
import shutil

#
#
#
def find_closest_xy_positions(dataset, well_position=(0.0, 0.0), num_positions=5):
    coords = dataset['Coord'].values
    x_coords = coords[:, 0]
    y_coords = coords[:, 1]
    
    # Find unique (x,y) positions
    xy_positions = np.column_stack((x_coords, y_coords))
    unique_xy, _ = np.unique(xy_positions, axis=0, return_inverse=True)
    
    # Calculate distances from well
    well_x, well_y = well_position
    distances = np.sqrt((unique_xy[:, 0] - well_x)**2 + (unique_xy[:, 1] - well_y)**2)
    
    # Get closest positions
    closest_indices = np.argsort(distances)[:num_positions]
    closest_positions = []
    
    # Generate colors and intensities (stronger for closer positions)
    color_intensities = np.linspace(1.0, 0.3, num_positions)
    base_colors = plt.cm.plasma(np.linspace(0.1, 0.9, num_positions))
    
    for i, idx in enumerate(closest_indices):
        x_pos, y_pos = unique_xy[idx]
        distance = distances[idx]
        
        closest_positions.append({
            'x_pos': x_pos,
            'y_pos': y_pos,
            'distance': distance,
            'color': 'k', #base_colors[i],
            'intensity': color_intensities[i],
            'label': f'({x_pos:.1f},{y_pos:.1f}) d={distance:.1f}m'
        })
    
    return closest_positions

#
#
#
def extract_stress_data_at_time(dataset, time_idx):
    coords = dataset['Coord'].values
    sxx_magnitude = dataset['Total Stress'][time_idx, :, 0].values / 1e6
    total_stress = dataset['Total Stress'][time_idx, :, :].values
    
    return {
        'sxx_magnitude': sxx_magnitude,
        'sigmaxx': total_stress[:, 0] / 1e6,
        'sigmayy': total_stress[:, 4] / 1e6,
        'sigmazz': total_stress[:, 8] / 1e6,
        'coords': coords
    }

#
#
#
def extract_vertical_profile(stress_data, xy_position, xy_tolerance=1e-6):
    coords = stress_data['coords']
    x_coords = coords[:, 0]
    y_coords = coords[:, 1]
    z_coords = coords[:, 2]
    
    # Find points at this (x,y) position
    mask = ((np.abs(x_coords - xy_position['x_pos']) <= xy_tolerance) & 
            (np.abs(y_coords - xy_position['y_pos']) <= xy_tolerance))
    
    if not np.any(mask):
        return None
    
    # Extract and sort by z-coordinate
    z_slice = z_coords[mask]
    sxx_slice = stress_data['sxx_magnitude'][mask]
    sigmazz_slice = stress_data['sigmazz'][mask]
    
    sort_idx = np.argsort(z_slice)
    
    return {
        'z_coords': z_slice[sort_idx],
        'sxx_magnitude': sxx_slice[sort_idx],
        'sigmazz_magnitude': sigmazz_slice[sort_idx],
        'xy_position': xy_position
    }


#
#
#
def extract_time_series_at_position(dataset, xy_position, target_depths, xy_tolerance=1e-6, z_tolerance=2.0):
    coords = dataset['Coord'].values
    x_coords = coords[:, 0]
    y_coords = coords[:, 1]
    z_coords = coords[:, 2]
    time_values = dataset.time.values
    
    time_series_data = {}
    
    for target_z in target_depths:
        # Find points at this position and depth
        mask = ((np.abs(x_coords - xy_position['x_pos']) <= xy_tolerance) & 
                (np.abs(y_coords - xy_position['y_pos']) <= xy_tolerance) & 
                (np.abs(z_coords - target_z) <= z_tolerance))
        
        if not np.any(mask):
            continue
        
        point_indices = np.where(mask)[0]
        
        # Extract time series
        sxx_series = []
        sigmaxx_series = []
        sigmayy_series = []
        sigmazz_series = []
        
        for t_idx in range(len(time_values)):
            total_stress = dataset['Total Stress'][t_idx, point_indices, :].values
            sxx_magnitude = dataset['Total Stress'][t_idx, point_indices, 0].values
            
            # Average over points at this position/depth
            sxx_series.append(np.nanmean(sxx_magnitude / 1e6))
            sigmaxx_series.append(np.nanmean(total_stress[:, 0] / 1e6))
            sigmayy_series.append(np.nanmean(total_stress[:, 4] / 1e6))
            sigmazz_series.append(np.nanmean(total_stress[:, 8] / 1e6))
        
        time_series_data[target_z] = {
            'sxx': np.array(sxx_series),
            'sigmaxx': np.array(sigmaxx_series),
            'sigmayy': np.array(sigmayy_series),
            'sigmazz': np.array(sigmazz_series)
        }
    
    return time_series_data

#
#
#
def find_closest_profile_for_depth(profiles, target_z):
    best_profile = None
    best_sxx = None
    best_sigmazz = None
    min_distance = float('inf')
    
    for profile in profiles:
        if profile is None:
            continue
            
        z_coords = profile['z_coords']
        z_distances = np.abs(z_coords - target_z)
        closest_idx = np.argmin(z_distances)
        closest_distance = z_distances[closest_idx]
        
        if closest_distance < min_distance:
            min_distance = closest_distance
            best_profile = profile
            
            # Get values at closest point or interpolate
            if closest_distance < 0.5:
                best_sxx = profile['sxx_magnitude'][closest_idx]
                best_sigmazz = profile['sigmazz_magnitude'][closest_idx]
            else:
                # Interpolate if we have enough points
                if len(z_coords) > 1:
                    best_sxx = np.interp(target_z, z_coords, profile['sxx_magnitude'])
                    best_sigmazz = np.interp(target_z, z_coords, profile['sigmazz_magnitude'])
                else:
                    best_sxx = profile['sxx_magnitude'][closest_idx]
                    best_sigmazz = profile['sigmazz_magnitude'][closest_idx]
    
    return best_profile, best_sxx, best_sigmazz

#
#
#
def calculate_stress_ranges(time_series_data):
    all_stress_values = []
    
    for depth_data in time_series_data.values():
        if isinstance(depth_data, dict):
            for key in ['sxx', 'sigmaxx', 'sigmayy', 'sigmazz']:
                if key in depth_data:
                    all_stress_values.extend(depth_data[key])
    
    if not all_stress_values:
        return (-100, -70)  # Default range
    
    stress_min = np.min(all_stress_values) - 7
    stress_max = np.max(all_stress_values) + 7
    
    return stress_min, stress_max

#
#
#
def create_profile_legend_elements(tracked_depths):
    legend_elements = []
    
    # Stress type elements
    legend_elements.append(Line2D([0], [0], color='black', linewidth=2, label=r'$\sigma_h$'))
    legend_elements.append(Line2D([0], [0], color='black', linewidth=2, linestyle='--', label=r'$\sigma_{zz}$'))
    
    # Depth elements
    for depth in tracked_depths:
        if depth > 0:
            label = f'Reservoir Z={depth:+.1f}m'
            color = 'red'
        else:
            label = f'Caprock Z={depth:+.1f}m'
            color = 'blue'
        legend_elements.append(Line2D([0], [0], color=color, linewidth=2, label=label))
    
    # Marker elements
    legend_elements.append(Line2D([0], [0], marker='o', color='black', linewidth=0, 
                                markersize=8, markerfacecolor='white', 
                                markeredgewidth=2, label=r'$\sigma_h$ markers'))
    legend_elements.append(Line2D([0], [0], marker='s', color='black', linewidth=0, 
                                markersize=8, markerfacecolor='white', 
                                markeredgewidth=2, label=r'$\sigma_{zz}$ markers'))
    
    return legend_elements

#
#
#
def calculate_stress_differences(time_series_data, tracked_depths, time_idx):
    if len(tracked_depths) != 2 or not all(depth in time_series_data for depth in tracked_depths):
        return None
    
    depth1, depth2 = tracked_depths  # +5m (reservoir), -5m (caprock)
    depth1_data = time_series_data[depth1]
    depth2_data = time_series_data[depth2]
    
    # Calculate Δσ = Stress - σzz at each depth
    sxx_reservoir = depth1_data['sxx'][time_idx]
    sigmazz_reservoir = depth1_data['sigmazz'][time_idx]
    delta_sigma_reservoir = sxx_reservoir - sigmazz_reservoir
    
    sxx_caprock = depth2_data['sxx'][time_idx]
    sigmazz_caprock = depth2_data['sigmazz'][time_idx]
    delta_sigma_caprock = sxx_caprock - sigmazz_caprock
    
    # Difference between caprock and reservoir Δσ
    delta_between = delta_sigma_caprock - delta_sigma_reservoir
    
    return {
        'delta_between': delta_between,
        'delta_reservoir': delta_sigma_reservoir,
        'delta_caprock': delta_sigma_caprock
    }


#
#
#
def create_video_with_ffmpeg(output_dir, output_filename, fps, cleanup_frames=True):
    print(f"Creating video with FFmpeg at {fps} fps...")

    # FFmpeg commands with GPU acceleration (try NVIDIA first, fallback to CPU)
    ffmpeg_commands = [
        # NVIDIA GPU encoding
        [
            'ffmpeg', '-y', '-framerate', str(fps),
            '-i', os.path.join(output_dir, 'frame_%04d.png'),
            '-c:v', 'h264_nvenc',
            '-preset', 'slow',
            '-crf', '10',
            '-pix_fmt', 'yuv420p',
            '-metadata', 'artist=Stress Combined Animation',
            output_filename
        ],
        # CPU fallback
        [
            'ffmpeg', '-y', '-framerate', str(fps),
            '-i', os.path.join(output_dir, 'frame_%04d.png'),
            '-c:v', 'libx264',
            '-preset', 'fast',
            '-crf', '18',
            '-pix_fmt', 'yuv420p',
            '-metadata', 'artist=Stress Combined Animation',
            output_filename
        ]
    ]

    # Try GPU encoding first, fallback to CPU
    success = False
    for i, cmd in enumerate(ffmpeg_commands):
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            encoder_type = "GPU (NVENC)" if i == 0 else "CPU"
            print(f"Video created successfully using {encoder_type} encoding!")
            success = True
            break
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            if i == 0:
                print("GPU encoding failed, trying CPU encoding...")
            else:
                print(f"FFmpeg failed: {e}")
                if hasattr(e, 'stderr') and e.stderr:
                    print(f"FFmpeg stderr: {e.stderr}")

    if not success:
        print("Both GPU and CPU encoding failed. Please check FFmpeg installation.")
        return False

    # Clean up temporary frames
#     if cleanup_frames:
#         try:
#             shutil.rmtree(output_dir)
#             print("Temporary frames cleaned up")
#         except Exception as e:
#             print(f"Warning: Could not clean up temporary frames: {e}")

    return True
