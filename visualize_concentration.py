import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from read_input import read_input
from read_input import read_vertex_file
from numba import njit, prange
import os


@njit(fastmath=True)
def green_function_2d(r_vec_sq, tau, pe):
    """
    Calculates the 2D Green's function for concentration.
    Args:
        r_vec_sq (float): Squared distance between target point and source point.
        tau (float): Time difference (t - t').
        pe (float): Peclet number.
    Returns:
        float: Concentration contribution from the source at this time step.
    """
    if tau <= 1e-12:
        return 0.0
    return (pe / (4 * np.pi * tau)) * np.exp(-r_vec_sq * pe / (4 * tau))


@njit(fastmath=True)
def green_function_3d(r_vec_sq, tau, pe, sigma=1.0):
    """
    Calculates the 3D Green's function for concentration.
    Args:
        r_vec_sq (float): Squared distance between target point and source point.
        tau (float): Time difference (t - t').
        pe (float): Peclet number.
    Returns:
        float: Concentration contribution from the source at this time step.
    """
    if tau <= 1e-12:
        return 0.0
    return (sigma * (pe / (4 * np.pi * tau))**(1.5)) * np.exp(-r_vec_sq * pe / (4 * tau))


@njit(parallel=True)
def calculate_concentration_field_2d(grid_points_x, grid_points_y, trajectory, t1, dt, pe):
    """
    Calculates the concentration for every point on a 2D grid.
    """
    concentration_field = np.zeros((len(grid_points_y), len(grid_points_x)))

    time_steps_to_consider = int(t1 / dt)
    if time_steps_to_consider > len(trajectory):
        time_steps_to_consider = len(trajectory)

    # Loop over each grid point in parallel
    for i in prange(len(grid_points_y)):
        for j in prange(len(grid_points_x)):
            target_point = np.array([grid_points_x[j], grid_points_y[i]])
            concentration = 0.0

            # Integrate over the particle's history up to time t1
            for step in range(time_steps_to_consider):
                source_pos = trajectory[step]
                tau = t1 - (step * dt)
                r_vec_sq = np.sum((target_point - source_pos) ** 2)

                concentration += green_function_2d(r_vec_sq, tau, pe) * dt

            concentration_field[i, j] = concentration

    return concentration_field


@njit(parallel=True)
def calculate_concentration_field_3d_slice(grid_points_x, grid_points_y, z_slice, trajectory, t1, dt, pe):
    """
    Calculates the concentration for every point on a 2D slice of a 3D grid.
    """
    concentration_field = np.zeros((len(grid_points_y), len(grid_points_x)))

    time_steps_to_consider = int(t1 / dt)
    if time_steps_to_consider > len(trajectory):
        time_steps_to_consider = len(trajectory)

    # Loop over each grid point in parallel
    for i in prange(len(grid_points_y)):
        for j in prange(len(grid_points_x)):
            target_point = np.array([grid_points_x[j], grid_points_y[i], z_slice])
            concentration = 0.0

            # Integrate over the particle's history up to time t1
            for step in range(time_steps_to_consider):
                # For 3D trajectory, we only need the position part (first 3 elements)
                source_pos = trajectory[step, :3]
                tau = t1 - (step * dt)
                r_vec_sq = np.sum((target_point - source_pos) ** 2)

                concentration += green_function_3d(r_vec_sq, tau, pe) * dt

            concentration_field[i, j] = concentration

    return concentration_field


@njit
def quaternion_to_rotation_matrix(q):
    """Converts a quaternion into a 3x3 rotation matrix."""
    w, x, y, z = q
    w2, x2, y2, z2 = w * w, x * x, y * y, z * z
    wx, wy, wz = w * x, w * y, w * z
    xy, xz, yz = x * y, x * z, y * z

    rot_matrix = np.empty((3, 3))
    rot_matrix[0, 0] = w2 + x2 - y2 - z2
    rot_matrix[0, 1] = 2 * (xy - wz)
    rot_matrix[0, 2] = 2 * (xz + wy)
    rot_matrix[1, 0] = 2 * (xy + wz)
    rot_matrix[1, 1] = w2 - x2 + y2 - z2
    rot_matrix[1, 2] = 2 * (yz - wx)
    rot_matrix[2, 0] = 2 * (xz - wy)
    rot_matrix[2, 1] = 2 * (yz + wx)
    rot_matrix[2, 2] = w2 - x2 - y2 + z2
    return rot_matrix


@njit(parallel=True, fastmath=True)
def _precompute_rotated_nodes(n_history_steps, n_nodes, trajectory_history, nodes_body_frame):
    """
    Precomputes the world coordinates of surface nodes for all history steps.
    """
    rotated_nodes_history = np.empty((n_history_steps, n_nodes, 3))
    for i in prange(n_history_steps):
        pos_center_hist = trajectory_history[i, :3]
        orient_hist = trajectory_history[i, 3:]
        rot_matrix_hist = quaternion_to_rotation_matrix(orient_hist)
        # Use matrix multiplication for all nodes at once for efficiency
        rotated_nodes_history[i, :, :] = pos_center_hist + (nodes_body_frame @ rot_matrix_hist.T)
    return rotated_nodes_history


# --- Concentration Field Calculation Functions ---

@njit(parallel=True)
def calculate_concentration_field_3d(grid_x, grid_y, grid_z, trajectory, t1, dt, pe):
    """Calculates concentration on a full 3D grid for a point particle."""
    conc_field = np.zeros((len(grid_x), len(grid_y), len(grid_z)))
    time_steps_to_consider = min(int(t1 / dt), len(trajectory))

    for i in prange(len(grid_x)):
        for j in prange(len(grid_y)):
            for k in prange(len(grid_z)):
                target_point = np.array([grid_x[i], grid_y[j], grid_z[k]])
                concentration = 0.0
                for step in range(time_steps_to_consider):
                    source_pos = trajectory[step, :3]
                    tau = t1 - (step * dt)
                    r_vec_sq = np.sum((target_point - source_pos)**2)
                    concentration += green_function_3d(r_vec_sq, tau, pe) * dt
                conc_field[i, j, k] = concentration
    return conc_field


@njit(parallel=True)
def calculate_concentration_field_3d_janus(grid_x, grid_y, grid_z, trajectory, nodes_body_frame, surface_sources, t1, dt, pe):
    """
    Calculates concentration on a full 3D grid for a Janus particle using pre-computation.
    """
    time_steps_to_consider = min(int(t1 / dt), len(trajectory))
    sigma_values = surface_sources
    n_nodes = len(nodes_body_frame)

    # --- Pre-computation Step ---
    rotated_nodes = _precompute_rotated_nodes(time_steps_to_consider, n_nodes, trajectory[:time_steps_to_consider],
                                              nodes_body_frame)
    # --- End Pre-computation ---

    conc_field = np.zeros((len(grid_x), len(grid_y), len(grid_z)))

    for i in prange(len(grid_x)):
        for j in prange(len(grid_y)):
            for k in prange(len(grid_z)):
                target_point = np.array([grid_x[i], grid_y[j], grid_z[k]])
                concentration = 0.0

                # Main loop now uses the pre-computed node positions
                for step in range(time_steps_to_consider):
                    tau = t1 - (step * dt)
                    for s in range(n_nodes):
                        world_pos = rotated_nodes[step, s, :]
                        sigma = sigma_values[s]
                        r_vec_sq = np.sum((target_point - world_pos) ** 2)
                        concentration += green_function_3d(r_vec_sq, tau, pe, sigma) * dt
                conc_field[i, j, k] = concentration
    return conc_field


# --- Plotting Functions ---

def plot_2d_concentration(grid_x, grid_y, concentration_field, trajectory, t1, dt, output_file, title):
    """Creates and saves a 2D heatmap of the concentration field."""
    print("Generating 2D plot...")
    fig, ax = plt.subplots(figsize=(10, 8))

    # Use pcolormesh for accurate grid plotting
    im = ax.pcolormesh(grid_x, grid_y, concentration_field, shading='gouraud', cmap='viridis')
    fig.colorbar(im, ax=ax, label='Chemical Concentration C')

    # Overlay trajectory
    time_steps_to_plot = min(int(t1 / dt), len(trajectory))
    traj_to_plot = trajectory[:time_steps_to_plot, :2]
    ax.plot(traj_to_plot[:, 0], traj_to_plot[:, 1], 'w-', lw=1.5, alpha=0.8, label='Trajectory')

    # Mark start and end points
    ax.plot(traj_to_plot[0, 0], traj_to_plot[0, 1], 'go', markersize=8, label='Start')
    ax.plot(traj_to_plot[-1, 0], traj_to_plot[-1, 1], 'ro', markersize=8, label=f'End (t={t1})')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title(title)
    ax.legend()
    ax.set_aspect('equal', 'box')
    plt.tight_layout()

    plt.savefig(output_file)
    print(f"Plot saved to {output_file}")


def plot_3d_concentration(concentration_field, grid_x, grid_y, grid_z, trajectory, t1, dt, output_file, title):
    """Creates and saves a 3D scatter plot of the concentration field."""
    print("Generating 3D plot...")
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Find points with significant concentration
    threshold = concentration_field.mean() + 0.1 * concentration_field.std()
    mask = concentration_field > threshold

    if not np.any(mask):
        print("Warning: No concentration values above the threshold. The plot may be empty.")
        # Fallback to a lower threshold
        threshold = np.percentile(concentration_field[concentration_field > 0], 75)
        mask = concentration_field > threshold
        if not np.any(mask):
            print("Fallback threshold also resulted in no points. Skipping scatter plot.")
            return

    # Get grid coordinates for plotting
    X, Y, Z = np.meshgrid(grid_x, grid_y, grid_z, indexing='ij')
    x_coords = X[mask]
    y_coords = Y[mask]
    z_coords = Z[mask]
    colors = concentration_field[mask]

    # Normalize colors
    norm = plt.Normalize(vmin=colors.min(), vmax=colors.max())

    # Scatter plot
    sc = ax.scatter(x_coords, y_coords, z_coords, c=colors, cmap='viridis', norm=norm, s=15, alpha=0.7)

    # Add color bar
    fig.colorbar(sc, ax=ax, shrink=0.5, aspect=10, label='Chemical Concentration C')

    # Overlay trajectory
    time_steps_to_plot = min(int(t1 / dt), len(trajectory))
    traj_to_plot = trajectory[:time_steps_to_plot, :3]
    ax.plot(traj_to_plot[:, 0], traj_to_plot[:, 1], traj_to_plot[:, 2], 'r-', lw=2.0, label='Trajectory')
    ax.scatter(traj_to_plot[0, 0], traj_to_plot[0, 1], traj_to_plot[0, 2], c='lime', s=50, label='Start')
    ax.scatter(traj_to_plot[-1, 0], traj_to_plot[-1, 1], traj_to_plot[-1, 2], c='magenta', s=50, label=f'End (t={t1})')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(title)
    ax.legend()

    # Improve layout and view
    ax.view_init(elev=20., azim=-65)
    plt.tight_layout()

    plt.savefig(output_file)
    print(f"Plot saved to {output_file}")


def plot_3d_sliced_view(concentration_field, grid_x, grid_y, grid_z, trajectory, t1, dt, output_file, title):
    """Creates a composite plot with two 2D slices and a 3D context view."""
    print("Generating 3D sliced view...")

    # --- 1. Data Preparation for Slicing ---
    time_steps_to_plot = min(int(t1 / dt), len(trajectory))
    traj_to_plot = trajectory[:time_steps_to_plot, :3]

    center_point = np.mean(traj_to_plot, axis=0)
    x_slice_idx = np.argmin(np.abs(grid_x - center_point[0]))
    y_slice_idx = np.argmin(np.abs(grid_y - center_point[1]))

    yz_slice_data = concentration_field[x_slice_idx, :, :].T
    xz_slice_data = concentration_field[:, y_slice_idx, :].T

    c_min = min(yz_slice_data.min(), xz_slice_data.min())
    c_max = max(yz_slice_data.max(), xz_slice_data.max())
    norm = Normalize(vmin=c_min, vmax=c_max)

    # --- 2. Create Figure and Subplots ---
    fig = plt.figure(figsize=(22, 7))
    ax1 = fig.add_subplot(1, 3, 1)
    ax2 = fig.add_subplot(1, 3, 2)
    ax3 = fig.add_subplot(1, 3, 3, projection='3d')
    fig.suptitle(title, fontsize=16)

    # --- 3. Plot YZ Slice (Left Panel) ---
    ax1.pcolormesh(grid_y, grid_z, yz_slice_data, cmap='viridis', norm=norm, shading='gouraud')
    ax1.plot(traj_to_plot[:, 1], traj_to_plot[:, 2], 'w-', lw=1.5, alpha=0.7)
    ax1.scatter(traj_to_plot[-1, 1], traj_to_plot[-1, 2], c='red', s=40, zorder=5)
    ax1.set_title(f'YZ Cross-section at X = {grid_x[x_slice_idx]:.2f}')
    ax1.set_xlabel('Y')
    ax1.set_ylabel('Z')
    ax1.set_aspect('equal', 'box')

    # --- 4. Plot XZ Slice (Middle Panel) ---
    im = ax2.pcolormesh(grid_x, grid_z, xz_slice_data, cmap='viridis', norm=norm, shading='gouraud')
    ax2.plot(traj_to_plot[:, 0], traj_to_plot[:, 2], 'w-', lw=1.5, alpha=0.7)
    ax2.scatter(traj_to_plot[-1, 0], traj_to_plot[-1, 2], c='red', s=40, zorder=5)
    ax2.set_title(f'XZ Cross-section at Y = {grid_y[y_slice_idx]:.2f}')
    ax2.set_xlabel('X')
    ax2.set_ylabel('Z')
    ax2.set_aspect('equal', 'box')

    fig.colorbar(im, ax=ax2, orientation='vertical', fraction=0.046, pad=0.04, label='Concentration C')

    # --- 5. Plot 3D Context View (Right Panel) with Gradient Trajectory ---
    points = traj_to_plot.reshape(-1, 1, 3)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    # Create a continuous norm to map from time to color
    cmap = plt.get_cmap('plasma')
    t_values = np.linspace(0, t1, len(traj_to_plot))
    line_norm = Normalize(vmin=t_values.min(), vmax=t_values.max())

    lc = Line3DCollection(segments, cmap=cmap, norm=line_norm)
    lc.set_array(t_values)
    lc.set_linewidth(2)
    ax3.add_collection(lc)

    # Add Start and End markers
    ax3.scatter(traj_to_plot[0, 0], traj_to_plot[0, 1], traj_to_plot[0, 2], c='lime', s=50, label='Start (t=0)',
                depthshade=False)
    ax3.scatter(traj_to_plot[-1, 0], traj_to_plot[-1, 1], traj_to_plot[-1, 2], c='red', s=50, label=f'End (t={t1:.1f})',
                depthshade=False)

    yy_plane, zz_plane = np.meshgrid(grid_y, grid_z)
    xx_plane_val = np.full_like(yy_plane, grid_x[x_slice_idx])
    ax3.plot_surface(xx_plane_val, yy_plane, zz_plane, color='red', alpha=0.15, rstride=5, cstride=5)

    xx_plane, zz_plane = np.meshgrid(grid_x, grid_z)
    yy_plane_val = np.full_like(xx_plane, grid_y[y_slice_idx])
    ax3.plot_surface(xx_plane, yy_plane_val, zz_plane, color='green', alpha=0.15, rstride=5, cstride=5)

    ax3.set_title('3D Trajectory and Slice Positions')
    ax3.set_xlabel('X')
    ax3.set_ylabel('Y')
    ax3.set_zlabel('Z')
    ax3.legend()

    max_range = np.array(
        [grid_x.max() - grid_x.min(), grid_y.max() - grid_y.min(), grid_z.max() - grid_z.min()]).max() / 2.0
    mid_x, mid_y, mid_z = center_point[0], center_point[1], center_point[2]
    ax3.set_xlim(mid_x - max_range, mid_x + max_range)
    ax3.set_ylim(mid_y - max_range, mid_y + max_range)
    ax3.set_zlim(mid_z - max_range, mid_z + max_range)
    ax3.view_init(elev=20., azim=-65)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(output_file)
    print(f"Sliced view plot saved to {output_file}")


# --- Main Execution ---


def main():
    parser = argparse.ArgumentParser(description="Visualize chemical concentration from a chemotaxi simulation.")
    parser.add_argument('--input-file', dest='input_file', type=str, required=True,
                        help='Path to the simulation input file (e.g., outputname.inputfile).')
    parser.add_argument('--trajectory-file', dest='trajectory_file', type=str, required=True,
                        help='Path to the trajectory file (e.g., outputname.config).')
    parser.add_argument('--time', dest='t1', type=float, required=True,
                        help='The specific time t1 at which to visualize the concentration field.')
    parser.add_argument('--output-file', dest='output_file', type=str, default='concentration_field.png',
                        help='Name of the output plot file.')
    parser.add_argument('--resolution', type=int, default=100,
                        help='Grid resolution for the visualization.')
    parser.add_argument('--padding', type=float, default=5.0,
                        help='Padding around the trajectory for the visualization grid.')
    parser.add_argument('--save-dat', dest='save_dat', action='store_true',
                        help='If specified, save the concentration field to a .dat file.')

    args = parser.parse_args()

    # 1. Read input files
    print("Reading input files...")
    sim_params = read_input.ReadInput(args.input_file)
    trajectory = np.loadtxt(args.trajectory_file)

    if args.t1 > sim_params.n_steps * sim_params.dt:
        print(f"Warning: Specified time t1 ({args.t1}) is greater than the total simulation time "
              f"({sim_params.n_steps * sim_params.dt}). Visualizing at the last time step.")
        args.t1 = sim_params.n_steps * sim_params.dt

    print(f"Simulation Parameters: Domain={sim_params.domain}, Pe={sim_params.peclet_number}, dt={sim_params.dt}")

    # 2. Setup the visualization grid
    pos_data = trajectory[:, :int(sim_params.domain.strip('D'))]  # 2 for 2D, 3 for 3D

    min_coords = pos_data.min(axis=0) - args.padding
    max_coords = pos_data.max(axis=0) + args.padding

    grid_points_x = np.linspace(min_coords[0], max_coords[0], args.resolution)
    grid_points_y = np.linspace(min_coords[1], max_coords[1], args.resolution)

    # 3. Calculate the concentration field
    print(f"Calculating concentration field at t = {args.t1}...")
    if sim_params.domain == '2D':
        concentration_field = calculate_concentration_field_2d(grid_points_x, grid_points_y, trajectory,
                                                               args.t1, sim_params.dt, sim_params.peclet_number)
        print("Calculation complete.")
        if args.save_dat:
            dat_filename = os.path.splitext(args.output_file)[0] + '.dat'
            print(f"Saving concentration data to {dat_filename}...")
            np.savetxt(dat_filename, concentration_field, fmt='%.8e')
            print("Data saved.")

        title = f'Chemical Concentration at t={args.t1}\nPe={sim_params.peclet_number}, ' \
                f'$\\Lambda$={sim_params.mobility_alpha}'
        plot_2d_concentration(grid_points_x, grid_points_y, concentration_field, trajectory,
                              args.t1, sim_params.dt, args.output_file,title)

    elif sim_params.domain == '3D':
        print(f"Using 3D resolution: {args.resolution}x{args.resolution}x{args.resolution}. This may take a while.")
        grid_points_z = np.linspace(min_coords[2], max_coords[2], args.resolution)

        if getattr(sim_params, 'particle_type', 'non-janus') == 'janus':
            print("Janus particle detected. Loading surface distribution...")
            try:
                structure_file_name = sim_params.structure
                structure_ref_config = read_vertex_file.read_vertex_file(structure_file_name[0])
                surface_file = sim_params.chemical_distribution_file
                surface_sources = np.loadtxt(surface_file[0])
                print(f"Loaded {len(surface_sources)} sources from {surface_file}")

                # --- Normalization Start ---
                sigma_values = np.asarray(surface_sources, dtype=float)
                total_sigma = np.sum(sigma_values)

                if total_sigma > 1e-12:
                    surface_sources = sigma_values / total_sigma
                    print(
                        f"Normalized {len(surface_sources)} sources. Total sigma before: {total_sigma:.4f}, after: 1.0")
                else:
                    surface_sources = len(structure_ref_config)  # Set all to 0 if sum is 0
                    print("Warning: Total sigma of sources is zero. All emissions will be zero.")
                # --- End Normalization ---

            except (AttributeError, FileNotFoundError):
                print(
                    f"Error: particle_type is 'janus' but file is missing or not found in {args.input_file}.")
                return

            concentration_field = calculate_concentration_field_3d_janus(grid_points_x, grid_points_y, grid_points_z,
                                                                         trajectory, structure_ref_config,
                                                                         surface_sources, args.t1,
                                                                         sim_params.dt, sim_params.peclet_number)
        else:
            print("Non-Janus (point) particle detected.")
            concentration_field = calculate_concentration_field_3d(grid_points_x, grid_points_y, grid_points_z,
                                                                   trajectory, args.t1, sim_params.dt,
                                                                   sim_params.peclet_number)
        print("Calculation complete.")
        if args.save_dat:
            dat_filename = os.path.splitext(args.output_file)[0] + '.dat'
            print(f"Saving 3D concentration data to {dat_filename}...")
            np.save(dat_filename, concentration_field)  # Use .npy for efficiency
            print("Data saved as a binary .npy file.")

            title = f'3D Chemical Concentration at t={args.t1}\nPe={sim_params.peclet_number}, ' \
                    f'type={getattr(sim_params, "particle_type", "non-janus")}'
            plot_3d_sliced_view(concentration_field, grid_points_x, grid_points_y, grid_points_z,
                                trajectory, args.t1, sim_params.dt, args.output_file, title)

    else:
        raise ValueError("Domain must be either '2D' or '3D'.")


if __name__ == '__main__':
    main()
