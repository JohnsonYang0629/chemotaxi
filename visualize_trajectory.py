import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MaxNLocator
import argparse
import os
import sys


def plot_trajectory(trajectory_data, output_filename, dt=1.0):
    """
    Plots a 2D or 3D particle trajectory based on the data dimension and saves it.
    The color of the trajectory corresponds to the simulation time, calculated using dt.
    """
    is_3d = trajectory_data.shape[1] == 3
    points = trajectory_data.reshape(-1, 1, trajectory_data.shape[1])
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    time_values = np.arange(len(trajectory_data)) * dt

    if is_3d:
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, projection='3d')

        lc = Line3DCollection(segments, cmap='viridis', norm=plt.Normalize(time_values.min(), time_values.max()))
        lc.set_array(time_values)
        lc.set_linewidth(2)
        line = ax.add_collection3d(lc)

        cbar = fig.colorbar(line, ax=ax, shrink=0.7)
        cbar.set_label('Time', rotation=270, labelpad=15)

        ax.plot([trajectory_data[0, 0]], [trajectory_data[0, 1]], [trajectory_data[0, 2]],
                'go', markersize=10, label='Start')
        ax.plot([trajectory_data[-1, 0]], [trajectory_data[-1, 1]], [trajectory_data[-1, 2]],
                'ro', markersize=10, label='End')

        ax.set_xlabel('X Coordinate')
        ax.set_ylabel('Y Coordinate')
        ax.set_zlabel('Z Coordinate')
        ax.set_title('3D Particle Trajectory')
        ax.legend()
        ax.grid(True)

        # --- 设置坐标轴等比例和边距 ---
        x_min, x_max = trajectory_data[:, 0].min(), trajectory_data[:, 0].max()
        y_min, y_max = trajectory_data[:, 1].min(), trajectory_data[:, 1].max()
        z_min, z_max = trajectory_data[:, 2].min(), trajectory_data[:, 2].max()

        range_x = np.ptp(trajectory_data[:, 0])
        range_y = np.ptp(trajectory_data[:, 1])
        range_z = np.ptp(trajectory_data[:, 2])

        padding_x = range_x * 0.1
        padding_y = range_y * 0.1
        padding_z = range_z * 0.1

        if padding_x == 0: padding_x = 1
        if padding_y == 0: padding_y = 1
        if padding_z == 0: padding_z = 1

        ax.set_xlim(x_min - padding_x, x_max + padding_x)
        ax.set_ylim(y_min - padding_y, y_max + padding_y)
        ax.set_zlim(z_min - padding_z, z_max + padding_z)

        ax.set_box_aspect((range_x, range_y, range_z))
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4))

    else:
        # --- 2D 绘图 ---
        fig, ax = plt.subplots(figsize=(10, 8))
        lc = LineCollection(segments, cmap='viridis', norm=plt.Normalize(time_values.min(), time_values.max()))
        lc.set_array(time_values)
        lc.set_linewidth(2)
        line = ax.add_collection(lc)

        cbar = fig.colorbar(line, ax=ax)
        cbar.set_label('Time', rotation=270, labelpad=15)

        ax.plot(trajectory_data[0, 0], trajectory_data[0, 1], 'go', markersize=10, label='Start')
        ax.plot(trajectory_data[-1, 0], trajectory_data[-1, 1], 'ro', markersize=10, label='End')

        ax.set_xlabel('X Coordinate')
        ax.set_ylabel('Y Coordinate')
        ax.set_title('2D Particle Trajectory')
        ax.legend()
        ax.grid(True)
        ax.set_aspect('equal', adjustable='box')

    # --- Save the figure ---
    try:
        plt.tight_layout()
        plt.savefig(output_filename, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"Trajectory plot successfully saved to: {output_filename}")
    except Exception as e:
        print(f"Error saving the figure: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    """
    Main function: parses command-line arguments, loads data, and calls the plotting function.
    """
    parser = argparse.ArgumentParser(
        description="Plots a 2D or 3D trajectory from a file. Automatically detects the dimension.",
        epilog="Example: python plot_trajectory_cli.py ./data/run.txt --dt 0.01"
    )
    parser.add_argument("trajectory_file", help="Path to the trajectory data file (space-separated values per line)")
    parser.add_argument("--dt", type=float, default=1.0,
                        help="The time step value between data points for the color bar's time axis. Defaults to 1.0.")
    args = parser.parse_args()

    input_path = args.trajectory_file
    if not os.path.exists(input_path):
        print(f"Error: File not found -> {input_path}", file=sys.stderr)
        sys.exit(1)

    print(f"Loading file: {input_path}...")
    try:
        data = np.loadtxt(input_path)
        if data.ndim == 1:
            data = data.reshape(1, -1)
    except Exception as e:
        print(f"Error: Failed to load file. Ensure it is a plain text file with space-separated numbers.\nDetails: {e}",
              file=sys.stderr)
        sys.exit(1)

    num_columns = data.shape[1]
    if num_columns == 2:
        print("Detected 2D data (2 columns: x, y)")
        plot_data = data
    elif num_columns >= 3:
        print(f"Detected 3D data ({num_columns} columns)")
        print("    --> Using the first 3 columns (x, y, z) for plotting.")
        plot_data = data[:, :3]
    else:
        print(f"Error: Unsupported file format. The file has {num_columns} columns.", file=sys.stderr)
        print("   Please ensure the file has at least 2 columns.", file=sys.stderr)
        sys.exit(1)

    base_name = os.path.splitext(input_path)[0]
    output_filename = base_name + ".trajectory.png"
    print(f"Using dt = {args.dt} for time calculation.")
    plot_trajectory(plot_data, output_filename, dt=args.dt)


if __name__ == "__main__":
    main()
