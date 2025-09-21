#!/usr/bin/env -S python -i

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter, FFMpegWriter
from matplotlib.lines import Line2D
from timestr import format_time_duration
from netcdf import read_netcdf

# Enable LaTeX rendering
# plt.rcParams['text.usetex'] = True
# plt.rcParams['font.family'] = 'serif'

plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.default'] = 'regular'
plt.rcParams['font.size'] = 14           # Base font size
plt.rcParams['axes.labelsize'] = 16      # X and Y axis labels
plt.rcParams['axes.titlesize'] = 18      # Subplot titles
plt.rcParams['xtick.labelsize'] = 14     # X-axis tick labels
plt.rcParams['ytick.labelsize'] = 14     # Y-axis tick labels
plt.rcParams['legend.fontsize'] = 14     # Legend text
plt.rcParams['figure.titlesize'] = 20    # Main figure title

# Extract vertical stress profiles at evenly spaced Y positions
def extract_vertical_profiles(y_coords, z_coords, s3_magnitude, sigmazz_magnitude, num_lines=10):
    y_min, y_max = y_coords.min(), y_coords.max()
    y_line_positions = np.linspace(y_min, y_max, num_lines)
    profiles = []

    for i, y_pos in enumerate(y_line_positions):
        y_tolerance = (y_max - y_min) / (num_lines * 1.5)
        mask = np.abs(y_coords - y_pos) <= y_tolerance

        if np.any(mask):
            z_slice = z_coords[mask]
            s3_mag_slice = s3_magnitude[mask]
            sigmazz_mag_slice = sigmazz_magnitude[mask]
            sort_idx = np.argsort(z_slice)
            z_sorted = z_slice[sort_idx]
            s3_mag_sorted = s3_mag_slice[sort_idx]
            sigmazz_mag_sorted = sigmazz_mag_slice[sort_idx]

            profiles.append({
                'y_position': y_pos,
                'z_coords': z_sorted,
                's3_magnitude': s3_mag_sorted,
                'sigmazz_magnitude': sigmazz_mag_sorted,
                'line_label': f'Y={y_pos:.1f}'
            })

    return profiles

# Create combined animation with S3 profiles and time series
def create_combined_animation(filename='plane_yz.cd', max_timesteps=None, num_lines=8):
    print(f"Creating combined S3 profiles and time series animation from {filename}")
    print(f"Number of vertical profiles: {num_lines}")

    # Load dataset
    filepath = f"run/cdf/{filename}"
    print(f"Loading {filepath}...")
    ds = read_netcdf(filepath)

    # Get time values and coordinates
    time_values = ds.time.values
    coords = ds['Coord'].values
    y_coords = coords[:, 1]
    z_coords = coords[:, 2]
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

    y_min, y_max = y_coords.min(), y_coords.max()
    z_min, z_max = z_coords.min(), z_coords.max()
    print(f"Coordinate range: Y=[{y_min:.1f}, {y_max:.1f}], Z=[{z_min:.1f}, {z_max:.1f}]")

    # Pre-compute time series data for tracked depths
    print("Computing time series data for tracked depths...")
    tracked_depths = [5.0, -5.0]
    time_series_data = {}
    
    for target_z in tracked_depths:
        mask = np.abs(z_coords - target_z) <= 2.0
        if np.any(mask):
            point_indices = np.where(mask)[0]
            
            sigmaxx_series = []
            sigmayy_series = []
            sigmazz_series = []
            s3_series = []
            
            for t_idx in range(len(time_values)):
                total_stress = ds['Total Stress'][t_idx, point_indices, :].values
                s3_magnitude = ds['S3 Magnitude'][t_idx, point_indices].values
                
                sigmaxx_vals = total_stress[:, 0] / 1e6
                sigmayy_vals = total_stress[:, 4] / 1e6
                sigmazz_vals = total_stress[:, 8] / 1e6  # zz component (index 8)
                s3_vals = s3_magnitude / 1e6
                
                sigmaxx_series.append(np.nanmean(sigmaxx_vals))
                sigmayy_series.append(np.nanmean(sigmayy_vals))
                sigmazz_series.append(np.nanmean(sigmazz_vals))
                s3_series.append(np.nanmean(s3_vals))
            
            time_series_data[target_z] = {
                'sigmaxx': np.array(sigmaxx_series),
                'sigmayy': np.array(sigmayy_series),
                'sigmazz': np.array(sigmazz_series),
                's3': np.array(s3_series)
            }

    # Find S3 magnitude range for consistent plotting
    print("Analyzing S3 magnitude range...")
    all_s3_magnitudes = []
    for t_idx in [0, n_times//2, -1] if n_times > 2 else [0]:
        s3_magnitude = ds['S3 Magnitude'][time_indices[t_idx], :].values / 1e6
        valid_mag = s3_magnitude[~np.isnan(s3_magnitude)]
        all_s3_magnitudes.extend(valid_mag)

    s3_min = np.percentile(all_s3_magnitudes, 1)
    s3_max = np.percentile(all_s3_magnitudes, 99)
    padding_mpa = 5
    s3_min_padded = s3_min - padding_mpa
    s3_max_padded = s3_max + padding_mpa

    # Find stress ranges for shared axis
    all_stress_values = []
    for depth_data in time_series_data.values():
        all_stress_values.extend(depth_data['sigmaxx'])
        all_stress_values.extend(depth_data['sigmayy'])
        all_stress_values.extend(depth_data['sigmazz'])
        all_stress_values.extend(depth_data['s3'])
    
    stress_min = np.min(all_stress_values) - 2
    stress_max = np.max(all_stress_values) + 2

    # Use shared stress range for S3 profiles too
    s3_min_shared = stress_min
    s3_max_shared = stress_max

    print(f"S3 magnitude range: [{s3_min:.1f}, {s3_max:.1f} MPa]")
    print(f"S3 magnitude range (with 5 MPa padding): [{s3_min_padded:.1f}, {s3_max_padded:.1f} MPa]")

    # Process all timesteps for profiles
    print("Processing timesteps...")
    frame_data_list = []

    for frame_idx, time_idx in enumerate(time_indices):
        current_time = time_values_anim[frame_idx]
        s3_magnitude = ds['S3 Magnitude'][time_idx, :].values / 1e6
        
        # Extract sigmazz magnitude from Total Stress tensor (zz component, index 8)
        total_stress = ds['Total Stress'][time_idx, :, :].values
        sigmazz_magnitude = total_stress[:, 8] / 1e6  # Convert Pa to MPa
        
        profiles = extract_vertical_profiles(y_coords, z_coords, s3_magnitude, sigmazz_magnitude, num_lines)

        frame_data_list.append({
            'time_idx': time_idx,
            'time_value': current_time,
            'profiles': profiles
        })

        if frame_idx % 5 == 0:
            print(f"  Processed frame {frame_idx+1}/{n_times}")

    # Create figure with shared stress axis
    fig, (ax_profiles, ax_timeseries) = plt.subplots(2, 1, figsize=(10, 12))

    # Define colors - depth colors match between plots
    profile_colors = plt.cm.tab10(np.linspace(0, 1, num_lines))
    depth_colors = ['red', 'blue']  # Red=+5m (reservoir), Blue=-5m (caprock)
    stress_colors = ['purple', 'orange']  # S3, sigmazz

    # Animation update function
    def update_frame(frame_idx):
        frame_data = frame_data_list[frame_idx]
        current_time = frame_data['time_value']
        current_time_years = time_years_anim[frame_idx]
        time_idx = frame_data['time_idx']

        ax_profiles.clear()
        ax_timeseries.clear()

        # Add figure title on top left with current time
        fig.suptitle(f'Time: {format_time_duration(current_time)}', fontsize=16, x=0.02, y=0.98, ha='left', va='top')

        # TOP PLOT: S3 and sigmazz Profiles
        # Reference profiles (initial)
        frame0 = frame_data_list[0]
        for i, profile in enumerate(frame0['profiles']):
            z_coords_prof = profile['z_coords']
            s3_magnitude_prof = profile['s3_magnitude']
            ax_profiles.plot(s3_magnitude_prof, z_coords_prof, color='k', alpha=0.5, linewidth=2, ls='--')
            
            # Reference sigmazz profile
            sigmazz_magnitude_prof = profile['sigmazz_magnitude']
            ax_profiles.plot(sigmazz_magnitude_prof, z_coords_prof, color='k', alpha=0.5, linewidth=2, ls=':')

        # Current profiles
        for i, profile in enumerate(frame_data['profiles']):
            z_coords_prof = profile['z_coords']
            s3_magnitude_prof = profile['s3_magnitude']
            sigmazz_magnitude_prof = profile['sigmazz_magnitude']
            
            # Plot S3 profile
            ax_profiles.plot(s3_magnitude_prof, z_coords_prof, color='purple', linewidth=2, alpha=0.8)
            
            # Plot sigmazz profile
            ax_profiles.plot(sigmazz_magnitude_prof, z_coords_prof, color='orange', linewidth=2, alpha=0.8)

        # Add tracking markers at monitored depths
        for j, (target_z, depth_color) in enumerate(zip(tracked_depths, depth_colors)):
            if target_z in time_series_data:
                # S3 circles
                s3_current = time_series_data[target_z]['s3'][time_idx]
                ax_profiles.plot(s3_current, target_z, 'o', markersize=12, 
                               color=depth_color, markerfacecolor='white', 
                               markeredgewidth=3, markeredgecolor=depth_color)
                
                # Sigmazz squares
                sigmazz_current = time_series_data[target_z]['sigmazz'][time_idx]
                ax_profiles.plot(sigmazz_current, target_z, 's', markersize=10, 
                               color=depth_color, markerfacecolor='white', 
                               markeredgewidth=3, markeredgecolor=depth_color)

        ax_profiles.set_ylabel('Z (m)', fontsize=12)
        ax_profiles.set_xlim(s3_min_shared, s3_max_shared)
        ax_profiles.set_ylim(z_max, z_min)  # Inverted Z axis
        ax_profiles.grid(True, alpha=0.3)
        ax_profiles.tick_params(axis='x', labelbottom=False)  # Remove x-axis labels
        
        # Add shared legend to top plot (upper right)
        legend_elements = []
        legend_elements.append(Line2D([0], [0], color='purple', linewidth=2, label=r'$S_3$'))
        legend_elements.append(Line2D([0], [0], color='orange', linewidth=2, label=r'$\sigma_{zz}$'))
        
        # Add depth legend elements
        for depth, depth_color in zip(tracked_depths, depth_colors):
            if depth > 0:
                label = f'Reservoir Z={depth:+.1f}m'
            else:
                label = f'Caprock Z={depth:+.1f}m'
            linestyle = '-' if depth > 0 else '--'
            legend_elements.append(Line2D([0], [0], color=depth_color, linestyle=linestyle, 
                                        linewidth=2, label=label))
        
        # Add marker legend elements
        legend_elements.append(Line2D([0], [0], marker='o', color='black', linewidth=0, 
                                    markersize=8, markerfacecolor='white', 
                                    markeredgewidth=2, label=r'$S_3$ markers'))
        legend_elements.append(Line2D([0], [0], marker='s', color='black', linewidth=0, 
                                    markersize=8, markerfacecolor='white', 
                                    markeredgewidth=2, label=r'$\sigma_{zz}$ markers'))
        
        ax_profiles.legend(handles=legend_elements, fontsize=9, loc='upper right', ncol=1)

        # BOTTOM PLOT: Time Series - Time on Y axis, Stress on X axis (inverted time)
        stress_components = ['s3', 'sigmazz']
        stress_labels = ['S3', 'σzz']
        stress_colors = ['purple', 'orange']
        
        for j, (target_z, depth_color) in enumerate(zip(tracked_depths, depth_colors)):
            if target_z in time_series_data:
                depth_data = time_series_data[target_z]
                
                for i, (component, label, stress_color) in enumerate(zip(stress_components, stress_labels, stress_colors)):
                    linestyle = '-' if target_z > 0 else '--'  # Solid for +5m, dashed for -5m
                    ax_timeseries.plot(depth_data[component], time_years, 
                                     color=stress_color, linestyle=linestyle, linewidth=2, alpha=0.8)
                    
                    # Highlight current time point
                    current_value = depth_data[component][time_idx]
                    marker_shape = 'o' if component == 's3' else 's'
                    ax_timeseries.plot(current_value, current_time_years, marker_shape, 
                                     markersize=10, color=depth_color, 
                                     markerfacecolor='white', markeredgewidth=3,
                                     markeredgecolor=depth_color)

        ax_timeseries.set_ylabel('Time (years)', fontsize=12)
        ax_timeseries.set_xlabel(r'Stress (MPa)', fontsize=12)
        ax_timeseries.set_ylim(time_years[-1], time_years[0])  # Invert time axis
        ax_timeseries.set_xlim(s3_min_shared, s3_max_shared)
        ax_timeseries.grid(True, alpha=0.3)
        
        # Add horizontal line at current time
        ax_timeseries.axhline(y=current_time_years, color='gray', linestyle='-', alpha=0.8, linewidth=2)
        
        # Calculate and display Δσ = S3 - σzz at each depth
        if len(tracked_depths) == 2 and all(depth in time_series_data for depth in tracked_depths):
            depth1, depth2 = tracked_depths  # +5m (reservoir, red), -5m (caprock, blue)
            depth1_data = time_series_data[depth1]  # Reservoir (red circles/squares)
            depth2_data = time_series_data[depth2]  # Caprock (blue circles/squares)
            
            # Calculate Δσ = S3 - σzz at each depth (circle - square)
            s3_reservoir = depth1_data['s3'][time_idx]        # Reservoir S3 (red circle)
            sigmazz_reservoir = depth1_data['sigmazz'][time_idx]  # Reservoir σzz (red square)
            delta_sigma_reservoir = s3_reservoir - sigmazz_reservoir  # S3 - σzz at reservoir
            
            s3_caprock = depth2_data['s3'][time_idx]       # Caprock S3 (blue circle)
            sigmazz_caprock = depth2_data['sigmazz'][time_idx] # Caprock σzz (blue square)
            delta_sigma_caprock = s3_caprock - sigmazz_caprock   # S3 - σzz at caprock
            
            # Difference between caprock and reservoir Δσ
            delta_between = delta_sigma_caprock - delta_sigma_reservoir
            
            # Create 3-line text box with clear naming and stress component indicators
            line1 = r'$\Delta\sigma$ Cap-Res (hh): ' + f'{delta_between:+.1f} MPa'
            line2 = r'$\Delta\sigma$ Caprock (hv): ' + f'{delta_sigma_caprock:+.1f} MPa'
            line3 = r'$\Delta\sigma$ Reservoir (hv): ' + f'{delta_sigma_reservoir:+.1f} MPa'
            
            diff_text = f'{line1}\n{line2}\n{line3}'
            
            # Position text above the gray line using data coordinates
            text_x = (s3_min_shared + s3_max_shared) / 2  # Center of stress range
            dtmax = time_years[-1] - time_years[0];
            
            text_y = time_years[0] + dtmax * 0.5  # Higher above current time
            ax_timeseries.text(text_x, text_y, diff_text, 
                             fontsize=12, ha='center', va='center',
                             bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9))

        return [ax_profiles, ax_timeseries]

    # Create animation
    print("Creating combined animation...")
    anim = FuncAnimation(fig, update_frame, frames=n_times,
                        interval=500, blit=False, repeat=True)

    plt.subplots_adjust(left=0.12, right=0.95, bottom=0.08, top=0.95, hspace=0.05)

    filename = f'combined_S3_animation.mp4'
    print(f"Saving {filename}...")

    writer = FFMpegWriter( fps=1000//300, metadata={}, bitrate=8000,
        extra_args=['-vcodec', 'libx264', '-pix_fmt', 'yuv420p', '-crf', '18', '-preset', 'slow'] )
    anim.save(filename, writer=writer, dpi=150)
    print("Animation saved!")

    return fig, anim

#
#
#
if __name__ == "__main__":
    filename = 'plane_yz.cd'

    # Create combined animation
    print("Creating combined S3 animation...")
    fig, anim = create_combined_animation(
        filename='plane_yz.cd',
        max_timesteps=10000,
        num_lines=1,
    )
