# visualize_trajectory.py

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.legend_handler import HandlerTuple
import argparse
import os
import sys
import scipy.io


def load_colormaps_from_mat_simplified(mat_file_path, indices=None):
    """
    Loads colormaps from a .mat file where they are stored as simple variables
    (e.g., 'color1', 'color2', ...), either all or a specific subset by index.
    """
    custom_colormaps = []
    try:
        mat_contents = scipy.io.loadmat(mat_file_path)

        if indices:
            print(f"Attempting to load specific colormap indices: {indices}")
            color_indices_to_load = indices
        else:
            num_colors = len([k for k in mat_contents if k.startswith('color')])
            print(f"Found {num_colors} colormaps in .mat file. Loading all of them.")
            color_indices_to_load = range(1, num_colors + 1)

        for i in color_indices_to_load:
            var_name = f'Colors{i}'
            if var_name in mat_contents:
                c_colors = mat_contents[var_name]
                if c_colors.max() > 1.0:
                    c_colors = c_colors / 255.0

                custom_cmap = LinearSegmentedColormap.from_list(f'custom_{var_name}', c_colors)
                custom_colormaps.append(custom_cmap)
                print(f"  -> Loaded '{var_name}'")
            else:
                print(f"  -> Warning: Colormap index {i} ('{var_name}') not found in .mat file.")

        if not custom_colormaps:
            print("Warning: No colormaps were loaded. Falling back to default colormaps.")
            return None
        return custom_colormaps

    except FileNotFoundError:
        print(f"Error: The colormap file was not found at '{mat_file_path}'", file=sys.stderr)
        return None
    except Exception as e:
        print(f"An error occurred while reading the .mat file: {e}", file=sys.stderr)
        return None


def load_trajectory_new_format(file_path):
    """
    Loads trajectory data for multiple particles from the block-based file format.
    """
    trajectories = {}
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
            line_idx = 0
            while line_idx < len(lines):
                try:
                    num_particles = int(lines[line_idx].strip())
                    line_idx += 1
                except (ValueError, IndexError):
                    line_idx += 1
                    continue
                for i in range(num_particles):
                    if line_idx < len(lines):
                        coords = np.array([float(p) for p in lines[line_idx].strip().split()])
                        particle_id = i + 1
                        if particle_id not in trajectories:
                            trajectories[particle_id] = []
                        trajectories[particle_id].append(coords)
                        line_idx += 1
                    else:
                        break
        for pid in trajectories:
            trajectories[pid] = np.array(trajectories[pid])
        return trajectories
    except FileNotFoundError:
        print(f"Error: Trajectory file not found at '{file_path}'", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred while reading the file: {e}", file=sys.stderr)
        sys.exit(1)


def plot_trajectories(trajectories, output_filename, dt=1.0, colormaps=None, center=None, plot_range=None,
                      solid_color=False):
    """
    Plots multiple 2D or 3D trajectories with options for solid or gradient colors.
    The legend displays the start and end markers for each particle.
    """
    if not trajectories:
        print("No trajectory data to plot.")
        return

    num_particles = len(trajectories)
    first_pid = next(iter(trajectories))
    is_3d = trajectories[first_pid].shape[1] > 4
    num_steps = len(trajectories[first_pid])
    time_values = np.arange(num_steps) * dt

    fig = plt.figure(figsize=(14, 10))
    ax_rect = [0.05, 0.05, 0.7, 0.9] if not solid_color and num_particles > 0 else [0.05, 0.05, 0.9, 0.9]

    if is_3d:
        ax = fig.add_axes(ax_rect, projection='3d')
        ax.set_title(f'{num_particles} Particle Trajectories (3D)')
        ax.set_xlabel('X Coordinate')
        ax.set_ylabel('Y Coordinate')
        ax.set_zlabel('Z Coordinate')
    else:
        ax = fig.add_axes(ax_rect)
        ax.set_title(f'{num_particles} Particle Trajectories (2D)')
        ax.set_xlabel('X Coordinate')
        ax.set_ylabel('Y Coordinate')
        ax.set_aspect('equal', 'box')

    if not colormaps:
        colormaps = [plt.get_cmap(name) for name in ['viridis', 'plasma', 'inferno', 'magma', 'cividis']]

    legend_handles = []

    for i, (pid, trajectory) in enumerate(trajectories.items()):
        cmap = colormaps[i % len(colormaps)]
        plot_data = trajectory[:, :3] if is_3d else trajectory[:, :2]

        if solid_color:
            line_color = cmap(0.65)
            if is_3d:
                ax.plot(plot_data[:, 0], plot_data[:, 1], plot_data[:, 2],
                        color=line_color, linewidth=2.5, alpha=0.8)
            else:
                ax.plot(plot_data[:, 0], plot_data[:, 1],
                        color=line_color, linewidth=2.5, alpha=0.8)
        else:
            points = plot_data.reshape(-1, 1, plot_data.shape[1])
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            norm = plt.Normalize(time_values.min(), time_values.max())

            if is_3d:
                lc = Line3DCollection(segments, cmap=cmap, norm=norm)
                lc.set_array(time_values)
                lc.set_linewidth(2.5)
                ax.add_collection3d(lc)
            else:
                lc = LineCollection(segments, cmap=cmap, norm=norm)
                lc.set_array(time_values)
                lc.set_linewidth(2.5)
                ax.add_collection(lc)

            cax = fig.add_axes([0.8, 0.1 + i * (0.8 / num_particles), 0.03, 0.8 / num_particles - 0.05])
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            fig.colorbar(sm, cax=cax, orientation='vertical', label=f'Time - Particle {pid}')

        start_color, end_color = cmap(0.0), cmap(1.0)

        # Plot markers without labels
        if is_3d:
            p_start = ax.scatter(plot_data[0, 0], plot_data[0, 1], plot_data[0, 2],
                                 color=start_color, s=120, marker='s', edgecolors='black', depthshade=False, zorder=5)
            p_end = ax.scatter(plot_data[-1, 0], plot_data[-1, 1], plot_data[-1, 2],
                               color=end_color, s=120, marker='o', edgecolors='black', depthshade=False, zorder=5)
        else:
            p_start = ax.scatter(plot_data[0, 0], plot_data[0, 1],
                                 color=start_color, s=120, marker='s', edgecolors='black', zorder=5)
            p_end = ax.scatter(plot_data[-1, 0], plot_data[-1, 1],
                               color=end_color, s=120, marker='o', edgecolors='black', zorder=5)

        # Create handles for the legend
        legend_handles.append(p_start)
        legend_handles.append(p_end)

    # Create labels for the legend
    legend_labels = []
    for i in range(1, num_particles + 1):
        legend_labels.append(f'Particle {i} Start')
        legend_labels.append(f'Particle {i} End')

    if center is not None and plot_range is not None:
        print("Using specified center and range for plot boundaries.")
        x0, y0, z0 = center
        Lx, Ly, Lz = plot_range
        ax.set_xlim(x0 - Lx / 2, x0 + Lx / 2)
        ax.set_ylim(y0 - Ly / 2, y0 + Ly / 2)
        if is_3d:
            ax.set_zlim(z0 - Lz / 2, z0 + Lz / 2)
            ax.set_box_aspect((Lx, Ly, Lz))
    else:
        print("Using automatic padding for plot boundaries.")
        all_coords = np.vstack(list(trajectories.values()))
        x_min, x_max = all_coords[:, 0].min(), all_coords[:, 0].max()
        y_min, y_max = all_coords[:, 1].min(), all_coords[:, 1].max()
        padding_x = (x_max - x_min) * 0.1 if (x_max - x_min) > 0 else 1
        padding_y = (y_max - y_min) * 0.1 if (y_max - y_min) > 0 else 1
        ax.set_xlim(x_min - padding_x, x_max + padding_y)
        ax.set_ylim(y_min - padding_y, y_max + padding_y)
        if is_3d:
            z_min, z_max = all_coords[:, 2].min(), all_coords[:, 2].max()
            padding_z = (z_max - z_min) * 0.1 if (z_max - z_min) > 0 else 1
            ax.set_zlim(z_min - padding_z, z_max + padding_z)
            ax.set_box_aspect((np.ptp(all_coords[:, 0]), np.ptp(all_coords[:, 1]), np.ptp(all_coords[:, 2])))

    ax.legend(legend_handles, legend_labels)
    ax.grid(False)

    try:
        plt.savefig(output_filename, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"Multi-trajectory plot successfully saved to: {output_filename}")
    except Exception as e:
        print(f"Error saving the figure: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    """
    Main function: parses arguments, loads data, and plots.
    """
    parser = argparse.ArgumentParser(
        description="Plots multiple 2D or 3D trajectories from a block-formatted simulation output file.",
        epilog="Examples:\n"
               "1. Default gradient plot: python visualize_trajectory.py results/run.config\n"
               "2. Solid color plot: python visualize_trajectory.py results/run.config --solid_color\n"
               "3. Select specific colormaps (e.g., 2, 4, 8): python visualize_trajectory.py results/run.config --cmap_indices 2 4 8\n"
               "4. Manual box: python visualize_trajectory.py results/run.config --center 0 0 5 --range 20 20 10",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("trajectory_file", help="Path to the trajectory data file (e.g., *.config)")
    parser.add_argument("--dt", type=float, default=1.0,
                        help="The time step value between data points. Defaults to 1.0.")
    parser.add_argument("--colormap_file", type=str, default="colormap.mat",
                        help="Path to the .mat file containing colormaps. Defaults to 'colormap.mat'.")
    parser.add_argument("--center", type=float, nargs=3, metavar=('X0', 'Y0', 'Z0'),
                        help="Specify the center of the plot box (x0 y0 z0). Requires --range.")
    parser.add_argument("--range", type=float, nargs=3, metavar=('Lx', 'Ly', 'Lz'),
                        help="Specify the dimensions (total length) of the plot box (Lx Ly Lz). Requires --center.")
    parser.add_argument("--solid_color", action='store_true',
                        help="Use a solid color for each trajectory instead of a time gradient. No colorbars will be drawn.")
    parser.add_argument("--cmap_indices", type=int, nargs='+',
                        help="A space-separated list of indices (e.g., 2 4 8) to select specific colormaps from the .mat file.")
    args = parser.parse_args()

    if (args.center and not args.range) or (args.range and not args.center):
        parser.error("--center and --range must be used together.")

    print(f"Attempting to load colormaps from: {args.colormap_file}")
    custom_colormaps = load_colormaps_from_mat_simplified(args.colormap_file, args.cmap_indices)

    print(f"Loading multi-particle trajectory file: {args.trajectory_file}...")
    trajectories = load_trajectory_new_format(args.trajectory_file)

    if not trajectories:
        print("No trajectories were loaded. Exiting.")
        sys.exit(1)

    print(f"Successfully loaded data for {len(trajectories)} particle(s).")

    base_name = os.path.splitext(args.trajectory_file)[0]
    output_filename = base_name + ".trajectories.png"

    plot_trajectories(trajectories, output_filename, dt=args.dt, colormaps=custom_colormaps,
                      center=args.center, plot_range=args.range, solid_color=args.solid_color)


if __name__ == "__main__":
    main()
