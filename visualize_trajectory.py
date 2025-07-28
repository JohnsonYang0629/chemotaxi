import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d import Axes3D
import argparse
import os
import sys


def plot_trajectory(trajectory_data, output_filename, dt=1.0):
    """
    Plots a 2D or 3D particle trajectory based on the data dimension and saves it.
    The color of the trajectory corresponds to the simulation time, calculated using dt.
    """

    # Check the dimension
    is_3d = trajectory_data.shape[1] == 3

    # Prepare segments for the gradient-colored line
    points = trajectory_data.reshape(-1, 1, trajectory_data.shape[1])
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    # Calculate the actual time values for color mapping
    time_values = np.arange(len(trajectory_data)) * dt

    if is_3d:
        # --- 3D Plotting ---
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, projection='3d')

        # For 3D, a scatter plot is used to simulate a gradient line for better clarity
        p = ax.scatter(trajectory_data[:, 0], trajectory_data[:, 1], trajectory_data[:, 2],
                       c=time_values, cmap='viridis', s=5)

        # Add a colorbar
        cbar = fig.colorbar(p, ax=ax, shrink=0.7)
        cbar.set_label('Time', rotation=270, labelpad=15)

        # Mark the start and end points
        ax.plot([trajectory_data[0, 0]], [trajectory_data[0, 1]], [trajectory_data[0, 2]],
                'go', markersize=10, label='Start')
        ax.plot([trajectory_data[-1, 0]], [trajectory_data[-1, 1]], [trajectory_data[-1, 2]],
                'ro', markersize=10, label='End')

        # Set plot properties
        ax.set_xlabel('X Coordinate')
        ax.set_ylabel('Y Coordinate')
        ax.set_zlabel('Z Coordinate')
        ax.set_title('3D Particle Trajectory')
        ax.legend()
        ax.grid(True)

    else:
        # --- 2D Plotting ---
        fig, ax = plt.subplots(figsize=(10, 8))

        # Use LineCollection to create a gradient-colored line
        lc = LineCollection(segments, cmap='viridis', norm=plt.Normalize(time_values.min(), time_values.max()))
        lc.set_array(time_values)
        lc.set_linewidth(2)
        line = ax.add_collection(lc)

        # Add a colorbar
        cbar = fig.colorbar(line, ax=ax)
        cbar.set_label('Time', rotation=270, labelpad=15)

        # Mark the start and end points
        ax.plot(trajectory_data[0, 0], trajectory_data[0, 1], 'go', markersize=10, label='Start')
        ax.plot(trajectory_data[-1, 0], trajectory_data[-1, 1], 'ro', markersize=10, label='End')

        # Set plot properties
        ax.set_xlabel('X Coordinate')
        ax.set_ylabel('Y Coordinate')
        ax.set_title('2D Particle Trajectory')
        ax.legend()
        ax.grid(True)
        ax.set_aspect('equal', adjustable='box')

    # --- Save the figure ---
    try:
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
    # 1. Set up the command-line argument parser
    parser = argparse.ArgumentParser(
        description="Plots a 2D or 3D trajectory from a file. Automatically detects the dimension.",
        epilog="Example: python plot_trajectory_cli.py ./data/run.txt --dt 0.01"
    )
    parser.add_argument("trajectory_file", help="Path to the trajectory data file (space-separated values per line)")
    parser.add_argument("--dt", type=float, default=1.0,
                        help="The time step value between data points for the color bar's time axis. Defaults to 1.0.")
    args = parser.parse_args()

    input_path = args.trajectory_file

    # 2. Check if the file exists
    if not os.path.exists(input_path):
        print(f"Error: File not found -> {input_path}", file=sys.stderr)
        sys.exit(1)

    print(f"Loading file: {input_path}...")

    # 3. Load the data
    try:
        data = np.loadtxt(input_path)
        if data.ndim == 1:  # If there is only one line of data, reshape it to a 2D array
            data = data.reshape(1, -1)
    except Exception as e:
        print(f"Error: Failed to load file. Ensure it is a plain text file with space-separated numbers.\nDetails: {e}",
              file=sys.stderr)
        sys.exit(1)

    num_columns = data.shape[1]

    # 4. Determine the dimension from the number of columns and prepare the data
    if num_columns == 2:
        print("Detected 2D data (2 columns: x, y)")
        plot_data = data
    elif num_columns == 7:
        print("Detected 3D data (7 columns: x, y, z, qw, qx, qy, qz)")
        print("    --> Using the first 3 columns (x, y, z) for plotting.")
        plot_data = data[:, :3]  # Extract only the first three columns
    else:
        print(f"Error: Unsupported file format. The file has {num_columns} columns.", file=sys.stderr)
        print("   Please ensure the 2D file has 2 columns and the 3D file has 7 columns.", file=sys.stderr)
        sys.exit(1)

    # 5. Generate the output filename and call the plotting function
    # Replace extensions like .txt, .dat with .trajectory.png
    base_name = os.path.splitext(input_path)[0]
    output_filename = base_name + ".trajectory.png"

    print(f"Using dt = {args.dt} for time calculation.")
    plot_trajectory(plot_data, output_filename, dt=args.dt)


if __name__ == "__main__":
    main()
