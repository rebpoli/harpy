#!/usr/bin/env -S python 

import os
import subprocess
import shutil
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.interpolate import griddata
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count
import time

from subsample import spatial_subsample
from netcdf import read_netcdf
from timestr import format_time_duration

def create_vector_segments(coords, vectors, pt_length=10):
    (center_x, center_y), (vec_x, vec_y) = coords[:, :2].T, vectors[:, :2].T
    return np.column_stack([center_x, center_y, vec_x, vec_y, np.full(len(center_x), pt_length)])

def subsample_data_indices(coords_2d, target_density):
    x_coords, y_coords = coords_2d[:, 0], coords_2d[:, 1]
    x_min, x_max = x_coords.min(), x_coords.max()
    y_min, y_max = y_coords.min(), y_coords.max()
    
    x_bins = np.linspace(x_min, x_max, target_density + 1)
    y_bins = np.linspace(y_min, y_max, target_density + 1)
    
    selected_indices = []
    for i in range(target_density):
        for j in range(target_density):
            x_mask = (x_coords >= x_bins[i]) & (x_coords < x_bins[i+1])
            y_mask = (y_coords >= y_bins[j]) & (y_coords < y_bins[j+1])
            candidates = np.where(x_mask & y_mask)[0]
            if len(candidates) > 0:
                selected_indices.append(candidates[0])
    
    return np.array(selected_indices)

def plot_segments(ax, segments, color, lengthscale, linewidth, alpha):
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    fig = ax.get_figure()
    fig_width_inch, fig_height_inch = fig.get_size_inches()
    
    bbox = ax.get_position()
    axis_width_points = bbox.width * fig_width_inch * 72
    axis_height_points = bbox.height * fig_height_inch * 72
    
    data_width, data_height = xlim[1] - xlim[0], ylim[1] - ylim[0]
    scale_x = data_width / axis_width_points
    scale_y = data_height / axis_height_points
    
    for center_x, center_y, vec_x, vec_y, pt_length in segments:
        if not (np.isnan(center_x) or np.isnan(center_y) or np.isnan(vec_x) or np.isnan(vec_y)):
            data_length_x = pt_length * scale_x
            data_length_y = pt_length * scale_y
            
            dx = vec_x * data_length_x / 2 * lengthscale
            dy = vec_y * data_length_y / 2 * lengthscale
            
            x_coords = [center_x - dx, center_x + dx]
            y_coords = [center_y - dy, center_y + dy]
            
            line = Line2D(x_coords, y_coords, color=color, linewidth=linewidth, alpha=alpha)
            ax.add_line(line)

def process_timestep_and_save_frame(args):
    """Single function that processes a timestep and saves the frame"""
    (time_idx, datasets, plot_configs, plot_setup_data, vmin, vmax,
     time_value, output_dir, dpi) = args
    
    # Process timestep data
    plot_data_list = []
    for config, setup_data in zip(plot_configs, plot_setup_data):
        ds = datasets[config['filename']]
        
        # Extract and subsample vectors
        s1_vectors = ds['S1'][time_idx, :, :].values[setup_data['subsample_indices']]
        s3_vectors = ds['S3'][time_idx, :, :].values[setup_data['subsample_indices']]
        
        x_comp_idx, y_comp_idx = setup_data['x_comp_idx'], setup_data['y_comp_idx']
        s1_2d = s1_vectors[:, [x_comp_idx, y_comp_idx]]
        s3_2d = s3_vectors[:, [x_comp_idx, y_comp_idx]]
        
        # Background field
        background_data = ds['VP Strain Rate'].sel(ten9_comp='zz')[time_idx, :].values
        
        # Interpolate background
        bg_xi = np.linspace(setup_data['x_min'], setup_data['x_max'], 100)
        bg_yi = np.linspace(setup_data['y_min'], setup_data['y_max'], 100)
        bg_Xi, bg_Yi = np.meshgrid(bg_xi, bg_yi)
        
        all_coords_2d = setup_data['coords'][:, [x_comp_idx, y_comp_idx]]
        bg_Zi = griddata(all_coords_2d, background_data, (bg_Xi, bg_Yi), method='linear')
        
        # Create vector segments
        subsampled_coords_2d = setup_data['coords'][setup_data['subsample_indices']][:, [x_comp_idx, y_comp_idx]]
        s1_segments = create_vector_segments(subsampled_coords_2d, s1_2d)
        s3_segments = create_vector_segments(subsampled_coords_2d, s3_2d)
        
        plot_data_list.append({
            'config': config,
            'setup_data': setup_data,
            'background': (bg_Xi, bg_Yi, bg_Zi),
            's1_segments': s1_segments,
            's3_segments': s3_segments,
        })
    
    # Create and save frame
    fig = plt.figure(figsize=(15, 5))
    gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 0.8, 0.08], wspace=0.2, left=0.03, right=0.95)
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    cbar_ax = fig.add_subplot(gs[0, 3])
    
    fig.suptitle(f'Time:{format_time_duration(time_value)}', fontsize=14, y=0.99)
    
    # Colorbar
    dummy_data = np.linspace(vmin, vmax, 100).reshape(10, 10)
    im = axes[0].imshow(dummy_data, cmap='RdYlBu_r', vmin=vmin, vmax=vmax)
    axes[0].clear()
    
    cbar = plt.colorbar(im, cax=cbar_ax)
    cbar.set_label('VP Strain Rate (1/s)', rotation=270, labelpad=15)
    cbar.ax.tick_params(labelsize=9)
    
    # Plot subplots
    for i, (ax, plot_data) in enumerate(zip(axes, plot_data_list)):
        ax.clear()
        
        # Background
        bg_Xi, bg_Yi, bg_Zi = plot_data['background']
        bg_Zi_clipped = np.clip(bg_Zi, vmin, vmax)
        ax.imshow(bg_Zi_clipped, extent=[bg_Xi.min(), bg_Xi.max(), bg_Yi.min(), bg_Yi.max()],
                 origin='lower', aspect='auto', interpolation='bilinear',
                 cmap='RdYlBu_r', vmin=vmin, vmax=vmax, alpha=0.6)
        
        # Vectors
        plot_segments(ax, plot_data['s3_segments'], 'lightgreen', 0.7, 1.5, 1)
        plot_segments(ax, plot_data['s1_segments'], 'black', 1.2, 1.5, 0.8)
        
        # Styling
        config, setup_data = plot_data['config'], plot_data['setup_data']
        ax.set_xlabel(config['x_axis'].upper())
        ax.set_ylabel(config['y_axis'].upper(), labelpad=-6)
        ax.set_title(config['title'])
        ax.set_xlim(setup_data['x_min'], setup_data['x_max'])
        ax.set_ylim(setup_data['y_min'], setup_data['y_max'])
        
        if i > 0:
            ax.invert_yaxis()
        
        # Legend
        legend_elements = [
            Line2D([0], [0], color='black', lw=2, label='S1'),
            Line2D([0], [0], color='lightgreen', lw=2, label='S3')
        ]
        ax.legend(handles=legend_elements, loc='upper right', fontsize=10)
    
    # Save frame
    frame_filename = os.path.join(output_dir, f'frame_{time_idx:04d}.png')
    fig.savefig(frame_filename, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    
    return frame_filename

def setup_plot_coordinates(datasets, plot_configs, vector_density):
    plot_setup_data = []
    for config in plot_configs:
        print(f"Setting up {config['filename']}...")
        
        ds = datasets[config['filename']]
        coords = ds['Coord'].values
        
        x_comp_idx = {'x': 0, 'y': 1, 'z': 2}[config['x_axis']]
        y_comp_idx = {'x': 0, 'y': 1, 'z': 2}[config['y_axis']]
        coords_2d = coords[:, [x_comp_idx, y_comp_idx]]
        
        x_min, x_max = coords_2d[:, 0].min(), coords_2d[:, 0].max()
        y_min, y_max = coords_2d[:, 1].min(), coords_2d[:, 1].max()
        
        subsample_indices = subsample_data_indices(coords_2d, vector_density)
        print(f"  Subsampled to: {len(subsample_indices)} points")
        
        plot_setup_data.append({
            'coords': coords, 'x_comp_idx': x_comp_idx, 'y_comp_idx': y_comp_idx,
            'subsample_indices': subsample_indices,
            'x_min': x_min, 'x_max': x_max, 'y_min': y_min, 'y_max': y_max
        })
    
    return plot_setup_data

def create_animated_vector_plots(max_timesteps=None, interval=300, vector_density=15, dpi=150):
    # Plot configurations
    plot_configs = [
        {'title': 'XY , 5m into the reservoir', 'filename': 'plane_xy.cd', 'x_axis': 'x', 'y_axis': 'y'},
        {'title': 'YZ , Well @ (0,0)', 'filename': 'plane_yz.cd', 'x_axis': 'y', 'y_axis': 'z'},
        {'title': 'YZ , Well @ (0,0)', 'filename': 'plane_yz_well.cd', 'x_axis': 'y', 'y_axis': 'z'}
    ]
    
    # Load datasets
    print("Loading datasets...")
    datasets = {}
    for config in plot_configs:
        filepath = f"run/cdf/{config['filename']}"
        ds = read_netcdf(filepath)
        datasets[config['filename']] = ds
        time_values = ds.time.values
        print(f"Found {len(time_values)} time steps")
    
    # Set color scale
    _data = ds['VP Strain Rate'].sel(ten9_comp='zz')
    vmin, vmax = float(_data.quantile(0.05).values), float(_data.quantile(0.95).values)
    
    n_times = min(max_timesteps, len(time_values))
    print(f"Animating {n_times} timesteps")
    print(f"Colorbar range: [{vmin:.3e}, {vmax:.3e}]")
    
    # Setup coordinates
    plot_setup_data = setup_plot_coordinates(datasets, plot_configs, vector_density)
    
    # Create output directory
    output_dir = 'temp_frames'
    os.makedirs(output_dir, exist_ok=True)
    
    # Single parallel processing: timestep processing + frame generation
    print(f"Processing timesteps and generating frames in parallel...")
    process_args = [
        (time_idx, datasets, plot_configs, plot_setup_data, vmin, vmax,
         time_values[time_idx], output_dir, dpi)
        for time_idx in range(n_times)
    ]
    
    max_workers = min(cpu_count() - 1, 12)
    completed = 0
    frame_files = [None] * n_times
    start_time = time.time()
    
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        future_to_idx = {
            executor.submit(process_timestep_and_save_frame, args): args[0] 
            for args in process_args
        }
        
        # Process completed tasks and show progress
        for future in as_completed(future_to_idx):
            time_idx = future_to_idx[future]
            try:
                frame_file = future.result()
                frame_files[time_idx] = frame_file
                completed += 1
                
                # Calculate progress
                progress_percent = (completed / n_times) * 100
                elapsed_time = time.time() - start_time
                
                print(f"\rProcessing: {completed}/{n_times} ({progress_percent:.1f}%) - using {max_workers} processors", end="", flush=True)
                
            except Exception as exc:
                print(f"\nFrame {time_idx} generated an exception: {exc}")
    
    print(f"\nGenerated {completed} frames in {time.time() - start_time:.1f}s")
    
    # Create video
    fps = 1000 // interval
    output_filename = 'S1_S3_maps.mp4'
    print(f"Creating video with FFmpeg at {fps} fps...")
    
    ffmpeg_commands = [
        # NVIDIA GPU encoding
        ['ffmpeg', '-y', '-framerate', str(fps), '-i', os.path.join(output_dir, 'frame_%04d.png'),
         '-c:v', 'h264_nvenc', '-preset', 'slow', '-crf', '18', '-pix_fmt', 'yuv420p',
         '-metadata', 'artist=S1/S3 Vector Animation', output_filename],
        # CPU fallback
        ['ffmpeg', '-y', '-framerate', str(fps), '-i', os.path.join(output_dir, 'frame_%04d.png'),
         '-c:v', 'libx264', '-preset', 'fast', '-crf', '18', '-pix_fmt', 'yuv420p',
         '-metadata', 'artist=S1/S3 Vector Animation', output_filename]
    ]
    
    success = False
    for i, cmd in enumerate(ffmpeg_commands):
        try:
            subprocess.run(cmd, capture_output=True, text=True, check=True)
            encoder_type = "GPU (NVENC)" if i == 0 else "CPU"
            print(f"Video created successfully using {encoder_type} encoding!")
            success = True
            break
        except (subprocess.CalledProcessError, FileNotFoundError):
            if i == 0:
                print("GPU encoding failed, trying CPU encoding...")
    
    if not success:
        print("Both GPU and CPU encoding failed. Please check FFmpeg installation.")
        return None
    
    # Cleanup
    try:
        shutil.rmtree(output_dir)
        print("Temporary frames cleaned up")
    except Exception as e:
        print(f"Warning: Could not clean up temporary frames: {e}")
    
    print(f"Animation saved as {output_filename}")
    return output_filename

if __name__ == "__main__":
    MAX_TIMESTEPS = 1500
    VECTOR_DENSITY = 15
    
    video_file = create_animated_vector_plots(
        max_timesteps=MAX_TIMESTEPS,
        vector_density=VECTOR_DENSITY,
        interval=200,
        dpi=150
    )
