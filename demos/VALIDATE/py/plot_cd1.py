#!/usr/bin/env -S python -i

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.animation import FuncAnimation, FFMpegWriter
from scipy.interpolate import griddata
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count

from subsample import spatial_subsample
from netcdf import read_netcdf
from timestr import format_time_duration

#
#
def create_vector_segments(coords, vectors, pt_length=10):
    (center_x, center_y), (vec_x, vec_y) = coords[:, :2].T, vectors[:, :2].T
    return  np.column_stack([ center_x, center_y, vec_x, vec_y, np.full(len(center_x), pt_length) ])

#
# Subsample points to get evenly distributed grid
#
def subsample_data_indices(coords_2d, target_density):
    x_coords = coords_2d[:, 0]
    y_coords = coords_2d[:, 1]

    x_min, x_max = x_coords.min(), x_coords.max()
    y_min, y_max = y_coords.min(), y_coords.max()

    # Create grid
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

#
# Process a single timestep for parallel execution
# 
def process_single_timestep(args):
    (time_idx, datasets, plot_configs, plot_setup_data, vmin, vmax,
     time_value) = args 

    frame_data = { 'time_idx': time_idx, 'time_value': time_value, 'plots': [] }

    for i, (config, setup_data) in enumerate(zip(plot_configs, plot_setup_data)):
        ds = datasets[config['filename']]

        # Extract data
        s1_vectors_full = ds['S1'][time_idx, :, :].values
        s3_vectors_full = ds['S3'][time_idx, :, :].values

        # Subsample vectors
        subsample_indices = setup_data['subsample_indices']
        s1_vectors = s1_vectors_full[subsample_indices]
        s3_vectors = s3_vectors_full[subsample_indices]

        # Extract 2D components for the plane
        x_comp_idx = setup_data['x_comp_idx']
        y_comp_idx = setup_data['y_comp_idx']
        s1_2d = s1_vectors[:, [x_comp_idx, y_comp_idx]]
        s3_2d = s3_vectors[:, [x_comp_idx, y_comp_idx]]

        # Get S3 Magnitude as background field (ONLY field we use)
        background_data = ds['S3 Magnitude'][time_idx, :].values / 1e6 # (MPa)

        # Create background interpolation
        bg_xi = np.linspace(setup_data['x_min'], setup_data['x_max'], 100)
        bg_yi = np.linspace(setup_data['y_min'], setup_data['y_max'], 100)
        bg_Xi, bg_Yi = np.meshgrid(bg_xi, bg_yi)

        coords = setup_data['coords']
        all_coords_2d = coords[:, [x_comp_idx, y_comp_idx]]
        bg_Zi = griddata(all_coords_2d, background_data, (bg_Xi, bg_Yi), method='linear')

        # Create vector segments with fixed screen size (10 points)
        subsampled_coords_2d = coords[subsample_indices][:, [x_comp_idx, y_comp_idx]]
        s1_segments = create_vector_segments(subsampled_coords_2d, s1_2d)
        s3_segments = create_vector_segments(subsampled_coords_2d, s3_2d)

        plot_data = {
            'config': config,
            'setup_data': setup_data,
            'background': (bg_Xi, bg_Yi, bg_Zi),
            's1_segments': s1_segments,
            's3_segments': s3_segments,
        }

        frame_data['plots'].append(plot_data)

    return frame_data

#
# Setup coordinate systems and subsampling for all plots
#
def setup_plot_coordinates(datasets, plot_configs, vector_density):
    plot_setup_data = []

    for config in plot_configs:
        print(f"Setting up {config['filename']}...")

        ds = datasets[config['filename']]
        coords = ds['Coord'].values

        # Component indices for the plane
        x_comp_idx = {'x': 0, 'y': 1, 'z': 2}[config['x_axis']]
        y_comp_idx = {'x': 0, 'y': 1, 'z': 2}[config['y_axis']]

        # Get 2D coordinates for this plane
        coords_2d = coords[:, [x_comp_idx, y_comp_idx]]

        # Find coordinate bounds
        x_min, x_max = coords_2d[:, 0].min(), coords_2d[:, 0].max()
        y_min, y_max = coords_2d[:, 1].min(), coords_2d[:, 1].max()

        print(f"  Range: X=[{x_min:.1f}, {x_max:.1f}], Y=[{y_min:.1f}, {y_max:.1f}]")
        print(f"  Data aspect ratio (Y/X): {(y_max-y_min)/(x_max-x_min):.3f}")

        # Subsample points
        subsample_indices = subsample_data_indices(coords_2d, vector_density)
        print(f"  Subsampled to: {len(subsample_indices)} points")

        setup_data = {
            'coords': coords,
            'x_comp_idx': x_comp_idx,
            'y_comp_idx': y_comp_idx,
            'subsample_indices': subsample_indices,
            'x_min': x_min, 'x_max': x_max,
            'y_min': y_min, 'y_max': y_max
        }

        plot_setup_data.append(setup_data)

    return plot_setup_data

#
# Plot line segments with fixed screen size using matplotlib lines
#
def plot_segments(ax, segments, color, linewidth, alpha):
    from matplotlib.lines import Line2D

    # Get current axis limits and figure size to calculate screen scaling
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    # Get figure size in inches and DPI
    fig = ax.get_figure()
    fig_width_inch, fig_height_inch = fig.get_size_inches()
    dpi = fig.dpi

    # Calculate axis size in points (there are 72 points per inch)
    bbox = ax.get_position()
    axis_width_points = bbox.width * fig_width_inch * 72
    axis_height_points = bbox.height * fig_height_inch * 72

    # Calculate data-to-screen scaling
    data_width = xlim[1] - xlim[0]
    data_height = ylim[1] - ylim[0]

    scale_x = data_width / axis_width_points
    scale_y = data_height / axis_height_points

    # Plot each segment
    for i in range(len(segments)):
        center_x, center_y, vec_x, vec_y, pt_length = segments[i]

        if not (np.isnan(center_x) or np.isnan(center_y) or np.isnan(vec_x) or np.isnan(vec_y)):
            # Convert screen length to data coordinates
            data_length_x = pt_length * scale_x
            data_length_y = pt_length * scale_y

            # Use the 2D vector components directly (already from unit 3D vector)
            # Scale by the screen length
            dx = vec_x * data_length_x / 2
            dy = vec_y * data_length_y / 2

            x_coords = [center_x - dx, center_x + dx]
            y_coords = [center_y - dy, center_y + dy]

            line = Line2D(x_coords, y_coords, color=color, linewidth=linewidth, alpha=alpha)
            ax.add_line(line)


#
#   Create animated S1 & S3 vector field plots with S3 Magnitude background
#
def create_animated_vector_plots(max_timesteps=None, interval=300,
                                vector_density=15):

    # Plot configurations
    plot_configs = [
        {'title': 'XY , 5m into the reservoir',    'filename': 'plane_xy.cd',         'x_axis': 'x',    'y_axis': 'y'},
        {'title': 'YZ , Well @ (0,0)',             'filename': 'plane_yz.cd',         'x_axis': 'y',    'y_axis': 'z'},
        {'title': 'YZ , Well @ (0,0)',             'filename': 'plane_yz_well.cd',    'x_axis': 'y',    'y_axis': 'z'}
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

    # Estimate reasonable vmin vmax
    vmin = float(ds['S3 Magnitude'].quantile(0.05).values)/1e6
    vmax = float(ds['S3 Magnitude'].quantile(0.95).values)/1e6

    # Determine timesteps to animate
    n_times = min(max_timesteps, len(time_values))
    print(f"Animating {n_times} timesteps")

    # Setup plot coordinates
    plot_setup_data = setup_plot_coordinates(datasets, plot_configs, vector_density)

    # Process all timesteps in parallel
    print("Processing timesteps...")
    process_args = []
    for time_idx in range(n_times):
        args = (time_idx, datasets, plot_configs, plot_setup_data, vmin, vmax,
                time_values[time_idx]) 
        process_args.append(args)
    with ProcessPoolExecutor(max_workers=min(cpu_count(), 6)) as executor:
        frame_data_list = list(executor.map(process_single_timestep, process_args))
    ##

    # Create figure with compact layout
    fig = plt.figure(figsize=(15, 5))
    gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 0.8, 0.08], wspace=0.2, left=0.03, right=0.95)
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    cbar_ax = fig.add_subplot(gs[0, 3])


    print("All plots set to square box aspect with minimal spacing")

    dummy_data = np.linspace(vmin, vmax, 100).reshape(10, 10)
    im = axes[0].imshow(dummy_data, cmap='RdYlBu_r', vmin=vmin, vmax=vmax)
    axes[0].clear()

    # Create the S3 Magnitude colorbar with proper spacing
    cbar = plt.colorbar(im, cax=cbar_ax)
    cbar.set_label('S3 (MPa)', rotation=270, labelpad=15)
    cbar.ax.tick_params(labelsize=9)

    print(f"Colorbar created with global S3 Magnitude range: [{vmin:.3e}, {vmax:.3e}]")

    # Animation update function
    def update_frame(frame_idx):
        frame_data = frame_data_list[frame_idx]
        current_time = frame_data['time_value']

        fig.suptitle(f'Time:{format_time_duration(current_time)}', fontsize=14, y=0.99)

        for i, (ax, plot_data) in enumerate(zip(axes, frame_data['plots'])):
            ax.clear()


            # S3 Magnitude background contour
            bg_Xi, bg_Yi, bg_Zi = plot_data['background']
            bg_Zi_clipped = np.clip(bg_Zi, vmin, vmax)

            # Use imshow for smooth linear interpolation
            im = ax.imshow(bg_Zi_clipped, extent=[bg_Xi.min(), bg_Xi.max(), bg_Yi.min(), bg_Yi.max()],
                          origin='lower', aspect='auto', interpolation='bilinear',
                          cmap='RdYlBu_r', vmin=vmin, vmax=vmax, alpha=0.6)

            # Plot S1 and S3 vectors with fixed screen size
            plot_segments(ax, plot_data['s3_segments'], 'green', 1.5, 0.7)
            plot_segments(ax, plot_data['s1_segments'], 'black', 1.5, 0.8)

            # Styling
            setup_data = plot_data['setup_data']
            config = plot_data['config']

#             if ( i < 2 ) : ax.set_box_aspect(1) 
            ax.set_xlabel(config['x_axis'].upper())
            ax.set_ylabel(config['y_axis'].upper(), labelpad=-6)  # Move Y label closer
            ax.set_title(config['title'])
            ax.set_xlim(setup_data['x_min'], setup_data['x_max'])
            ax.set_ylim(setup_data['y_min'], setup_data['y_max'])

            if i > 0 : ax.invert_yaxis();

            # Legend
            from matplotlib.lines import Line2D
            legend_elements = [
                Line2D([0], [0], color='black', lw=2, label='S1'),
                Line2D([0], [0], color='green', lw=2, label='S3')
            ]
            ax.legend(handles=legend_elements, loc='upper right', fontsize=10)

        return axes

    # Create animation
    print("Creating animation...")
    anim = FuncAnimation(fig, update_frame, frames=n_times,
                        interval=interval, blit=False, repeat=True)

    # Save animation as MP4 with high quality
    filename = f'S1_S3_maps.mp4'
    print(f"Saving {filename}...")

    # Use FFmpeg writer with high quality settings
    writer = FFMpegWriter(
        fps=1000//interval,           # Frame rate based on interval
        metadata=dict(artist='S1/S3 Vector Animation'),
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

    MAX_TIMESTEPS = 1500      # Number of timesteps to animate
    VECTOR_DENSITY = 15     # 15x15 grid of vectors

    fig, anim = create_animated_vector_plots(
        max_timesteps=MAX_TIMESTEPS,
        vector_density=VECTOR_DENSITY,
        interval=200  # 400ms between frames
    )

