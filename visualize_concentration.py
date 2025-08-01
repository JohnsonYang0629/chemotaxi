import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from read_input import read_input
from read_input import read_vertex_file
from numba import njit, prange
import os
import sys


# --- Calculation Functions ---
@njit(fastmath=True)
def green_function_2d(r_vec_sq, tau, pe):
    if tau <= 1e-12: return 0.0
    return (pe / (4 * np.pi * tau)) * np.exp(-r_vec_sq * pe / (4 * tau))


@njit(fastmath=True)
def green_function_3d(r_vec_sq, tau, pe, sigma=1.0):
    if tau <= 1e-12: return 0.0
    return (sigma * (pe / (4 * np.pi * tau)) ** (1.5)) * np.exp(-r_vec_sq * pe / (4 * tau))


@njit(parallel=True)
def calculate_concentration_field_2d(grid_points_x, grid_points_y, trajectory, t1, dt, pe):
    concentration_field = np.zeros((len(grid_points_y), len(grid_points_x)))
    time_steps_to_consider = min(int(t1 / dt), len(trajectory))
    for i in prange(len(grid_points_y)):
        for j in prange(len(grid_points_x)):
            target_point = np.array([grid_points_x[j], grid_points_y[i]])
            concentration = 0.0
            for step in range(time_steps_to_consider):
                source_pos = trajectory[step]
                tau = t1 - (step * dt)
                r_vec_sq = np.sum((target_point - source_pos) ** 2)
                concentration += green_function_2d(r_vec_sq, tau, pe) * dt
            concentration_field[i, j] = concentration
    return concentration_field


@njit
def quaternion_to_rotation_matrix(q):
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
    rotated_nodes_history = np.empty((n_history_steps, n_nodes, 3))
    for i in prange(n_history_steps):
        pos_center_hist = trajectory_history[i, :3]
        orient_hist = trajectory_history[i, 3:]
        rot_matrix_hist = quaternion_to_rotation_matrix(orient_hist)
        rotated_nodes_history[i, :, :] = pos_center_hist + (nodes_body_frame @ rot_matrix_hist.T)
    return rotated_nodes_history


@njit(parallel=True)
def calculate_concentration_field_3d(grid_x, grid_y, grid_z, trajectory, t1, dt, pe):
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
                    r_vec_sq = np.sum((target_point - source_pos) ** 2)
                    concentration += green_function_3d(r_vec_sq, tau, pe) * dt
                conc_field[i, j, k] = concentration
    return conc_field


@njit(parallel=True)
def calculate_concentration_field_3d_janus(grid_x, grid_y, grid_z, trajectory, nodes_body_frame, surface_sources, t1,
                                           dt, pe):
    time_steps_to_consider = min(int(t1 / dt), len(trajectory))
    n_nodes = len(nodes_body_frame)
    rotated_nodes = _precompute_rotated_nodes(time_steps_to_consider, n_nodes, trajectory[:time_steps_to_consider],
                                              nodes_body_frame)
    conc_field = np.zeros((len(grid_x), len(grid_y), len(grid_z)))
    for i in prange(len(grid_x)):
        for j in prange(len(grid_y)):
            for k in prange(len(grid_z)):
                target_point = np.array([grid_x[i], grid_y[j], grid_z[k]])
                concentration = 0.0
                for step in range(time_steps_to_consider):
                    tau = t1 - (step * dt)
                    for s in range(n_nodes):
                        world_pos = rotated_nodes[step, s, :]
                        sigma = surface_sources[s]
                        r_vec_sq = np.sum((target_point - world_pos) ** 2)
                        concentration += green_function_3d(r_vec_sq, tau, pe, sigma) * dt
                conc_field[i, j, k] = concentration
    return conc_field


# --- Plotting Functions ---
def plot_3d_sliced_view(concentration_field, grid_x, grid_y, grid_z, trajectory, t1, dt, output_file, title):
    """Creates a composite plot with two 2D slices and a 3D context view."""
    print("Generating 3D sliced view...")
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

    fig = plt.figure(figsize=(22, 7))
    ax1 = fig.add_subplot(1, 3, 1)
    ax2 = fig.add_subplot(1, 3, 2)
    ax3 = fig.add_subplot(1, 3, 3, projection='3d')
    fig.suptitle(title, fontsize=16)

    ax1.pcolormesh(grid_y, grid_z, yz_slice_data, cmap='viridis', norm=norm, shading='gouraud')
    ax1.plot(traj_to_plot[:, 1], traj_to_plot[:, 2], 'w-', lw=1.5, alpha=0.7)
    ax1.scatter(traj_to_plot[-1, 1], traj_to_plot[-1, 2], c='red', s=40, zorder=5)
    ax1.set_title(f'YZ Cross-section at X = {grid_x[x_slice_idx]:.2f}')
    ax1.set_xlabel('Y')
    ax1.set_ylabel('Z')
    ax1.set_aspect('equal', 'box')

    im = ax2.pcolormesh(grid_x, grid_z, xz_slice_data, cmap='viridis', norm=norm, shading='gouraud')
    ax2.plot(traj_to_plot[:, 0], traj_to_plot[:, 2], 'w-', lw=1.5, alpha=0.7)
    ax2.scatter(traj_to_plot[-1, 0], traj_to_plot[-1, 2], c='red', s=40, zorder=5)
    ax2.set_title(f'XZ Cross-section at Y = {grid_y[y_slice_idx]:.2f}')
    ax2.set_xlabel('X')
    ax2.set_ylabel('Z')
    ax2.set_aspect('equal', 'box')
    fig.colorbar(im, ax=ax2, orientation='vertical', fraction=0.046, pad=0.04, label='Concentration C')

    points = traj_to_plot.reshape(-1, 1, 3)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    t_values = np.linspace(0, t1, len(traj_to_plot))
    lc = Line3DCollection(segments, cmap=plt.get_cmap('plasma'),
                          norm=Normalize(vmin=t_values.min(), vmax=t_values.max()))
    lc.set_array(t_values)
    lc.set_linewidth(2)
    ax3.add_collection(lc)
    ax3.scatter(traj_to_plot[0, 0], traj_to_plot[0, 1], traj_to_plot[0, 2], c='lime', s=50, label='Start',
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
    plt.savefig(output_file, dpi=150)
    print(f"Sliced view plot saved to {output_file}")


# --- Main Execution Logic ---

def run_static_mode(args, sim_params, trajectory, grid_x, grid_y, grid_z):
    """Handles logic for generating a single static plot."""
    print(f"--- Running in STATIC mode for time t = {args.t1} ---")
    output_filename = args.output_file if args.output_file.lower().endswith('.png') else args.output_file + '.png'

    if sim_params.domain == '2D':
        print("Calculating 2D concentration field...")
        concentration_field = calculate_concentration_field_2d(grid_x, grid_y, trajectory,
                                                               args.t1, sim_params.dt, sim_params.peclet_number)
        print("Calculation complete.")
        title = f'Chemical Concentration at t={args.t1:.2f}\nPe={sim_params.peclet_number}, $\\Lambda$={sim_params.mobility_alpha}'

        fig, ax = plt.subplots(figsize=(10, 8))
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_title(title)
        ax.set_aspect('equal', 'box')

        im = ax.pcolormesh(grid_x, grid_y, concentration_field, shading='gouraud', cmap='viridis', zorder=1)
        ax.grid(True)
        fig.colorbar(im, ax=ax, label='Chemical Concentration C')

        time_steps_to_plot = min(int(args.t1 / sim_params.dt), len(trajectory))
        traj_to_plot = trajectory[:time_steps_to_plot]
        ax.plot(traj_to_plot[:, 0], traj_to_plot[:, 1], 'w-', lw=1.5, alpha=0.8, label='Trajectory', zorder=2)
        ax.plot(traj_to_plot[0, 0], traj_to_plot[0, 1], 'go', markersize=8, label='Start', zorder=3)
        ax.plot(traj_to_plot[-1, 0], traj_to_plot[-1, 1], 'ro', markersize=8, label=f'Current (t={args.t1:.2f})',
                zorder=3)
        ax.legend()
        plt.savefig(output_filename, dpi=150)
        print(f"Static 2D plot saved to {output_filename}")

    elif sim_params.domain == '3D':
        print(f"Calculating 3D concentration field for sliced view...")
        if getattr(sim_params, 'particle_type', 'non-janus') == 'janus':
            print("Janus particle detected. Loading surface distribution...")
            try:
                nodes = read_vertex_file.read_vertex_file(sim_params.structure[0])
                sources = np.loadtxt(sim_params.chemical_distribution_file[0])
                total_sigma = np.sum(sources)
                if total_sigma > 1e-12:
                    sources /= total_sigma
                print(f"Loaded and normalized {len(sources)} sources.")
                concentration_field = calculate_concentration_field_3d_janus(grid_x, grid_y, grid_z, trajectory, nodes,
                                                                             sources, args.t1, sim_params.dt,
                                                                             sim_params.peclet_number)
            except Exception as e:
                print(f"Error loading Janus particle data: {e}", file=sys.stderr)
                return
        else:
            print("Non-Janus (point) particle detected.")
            concentration_field = calculate_concentration_field_3d(grid_x, grid_y, grid_z, trajectory, args.t1,
                                                                   sim_params.dt, sim_params.peclet_number)

        print("Calculation complete.")
        title = f'3D Sliced View at t={args.t1:.2f}\n' \
                f'Pe={sim_params.peclet_number}, $\\Lambda$={sim_params.mobility_alpha}, ' \
                f'Type={getattr(sim_params, "particle_type", "non-janus")}'
        plot_3d_sliced_view(concentration_field, grid_x, grid_y, grid_z, trajectory, args.t1, sim_params.dt,
                            output_filename, title)


def run_animation_mode(args, sim_params, trajectory, grid_x, grid_y, grid_z):
    """Handles logic for generating an animation."""
    if sim_params.domain == '3D':
        print("Animation mode for 3D is not yet implemented. Please use static mode.", file=sys.stderr)
        return

    total_time = len(trajectory) * sim_params.dt
    print(f"--- Running in ANIMATION mode for total time = {total_time:.2f} ---")

    if args.writer == 'gif':
        output_filename = args.output_file if args.output_file.lower().endswith('.gif') else args.output_file + '.gif'
        writer_instance = animation.PillowWriter(fps=15)
        print("Using Pillow writer to generate GIF.")
    else:
        output_filename = args.output_file if args.output_file.lower().endswith('.mp4') else args.output_file + '.mp4'
        writer_instance = animation.FFMpegWriter(fps=15, metadata=dict(artist='Me'), bitrate=1800)
        print("Using FFMpeg writer to generate MP4.")

    animation_times = np.arange(0, total_time, args.frame_interval)
    if not animation_times.any() or animation_times[-1] < total_time:
        animation_times = np.append(animation_times, total_time)

    # --- FIX 1: Pre-calculate color limits to prevent flickering ---
    print("Pre-calculating color limits for stable animation...")
    final_field = calculate_concentration_field_2d(grid_x, grid_y, trajectory, total_time, sim_params.dt,
                                                   sim_params.peclet_number)
    vmin, vmax = final_field.min(), final_field.max()
    print(f"Global color range set to: [{vmin:.4f}, {vmax:.4f}]")

    # --- Setup initial plot ---
    fig, ax = plt.subplots(figsize=(10, 8))
    title_obj = ax.set_title('')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_aspect('equal', 'box')

    # Create the initial plot objects that will be updated each frame
    im = ax.pcolormesh(grid_x, grid_y, np.zeros_like(final_field), shading='gouraud', cmap='viridis', vmin=vmin,
                       vmax=vmax, zorder=1)
    line, = ax.plot([], [], 'w-', lw=1.5, alpha=0.8, zorder=2)
    start_dot, = ax.plot(trajectory[0, 0], trajectory[0, 1], 'go', markersize=8, label='Start', zorder=3)
    current_dot, = ax.plot([], [], 'ro', markersize=8, label='Current Position', zorder=3)
    ax.legend(loc='upper right')
    ax.grid(True)

    # --- FIX 2: Create the colorbar ONCE, outside the update loop ---
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label('Chemical Concentration C')

    def update(frame_time):
        # Calculate new data
        concentration_field = calculate_concentration_field_2d(grid_x, grid_y, trajectory, frame_time, sim_params.dt,
                                                               sim_params.peclet_number)

        # Update existing plot objects instead of clearing and redrawing
        im.set_array(concentration_field.ravel())

        current_step = min(int(frame_time / sim_params.dt), len(trajectory) - 1)
        traj_to_plot = trajectory[:current_step + 1]

        line.set_data(traj_to_plot[:, 0], traj_to_plot[:, 1])
        current_dot.set_data(traj_to_plot[-1, 0], traj_to_plot[-1, 1])

        title_obj.set_text(
            f't={frame_time:.2f} s; Pe={sim_params.peclet_number}, $\\Lambda$={sim_params.mobility_alpha}')

        print(f"Rendering frame for time: {frame_time:.2f} / {total_time:.2f}")
        return im, line, current_dot, title_obj

    print("Starting animation rendering... This may take a while.")
    try:
        ani = animation.FuncAnimation(fig, update, frames=animation_times, blit=False, repeat=False)
        ani.save(output_filename, writer=writer_instance)
        print(f"Animation successfully saved to {output_filename}")
    except FileNotFoundError:
        print(f"\n--- ERROR: '{args.writer}' writer not found ---", file=sys.stderr)
        if args.writer == 'ffmpeg':
            print("To save MP4 videos, please install ffmpeg.", file=sys.stderr)
        else:  # Pillow
            print("To save GIFs, please install the Pillow library: pip install Pillow", file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Visualize chemical concentration from a chemotaxi simulation statically or as an animation.")
    parser.add_argument('--input-file', dest='input_file', type=str, required=True,
                        help='Path to simulation input file.')
    parser.add_argument('--trajectory-file', dest='trajectory_file', type=str, required=True,
                        help='Path to trajectory file.')
    parser.add_argument('--output-file', dest='output_file', type=str, default='concentration_field',
                        help='Base name of output file (no extension).')
    parser.add_argument('--mode', choices=['static', 'animation'], default='static',
                        help="Choose visualization mode. Default: static.")
    parser.add_argument('--writer', choices=['ffmpeg', 'gif'], default='gif',
                        help="[ANIMATION mode] Writer to use for saving animation. 'ffmpeg' for MP4, 'gif' for GIF. Default: gif.")
    parser.add_argument('--time', dest='t1', type=float, help='[STATIC mode] Specific time t1 to visualize.')
    parser.add_argument('--frame-interval', type=float, default=5.0,
                        help='[ANIMATION mode] Simulation time between frames. Default: 5.0.')
    parser.add_argument('--resolution', type=int, default=100, help='Grid resolution for visualization.')
    parser.add_argument('--padding', type=float, default=5.0, help='Padding around the trajectory for the grid.')
    args = parser.parse_args()

    if args.mode == 'static' and args.t1 is None:
        parser.error("--time is required for --mode='static'")

    print("Reading input files...")
    sim_params = read_input.ReadInput(args.input_file)
    trajectory = np.loadtxt(args.trajectory_file)
    if trajectory.ndim == 1: trajectory = trajectory.reshape(1, -1)
    print(f"Simulation Parameters: Domain={sim_params.domain}, Pe={sim_params.peclet_number}, dt={sim_params.dt}")

    pos_data = trajectory[:, :2] if sim_params.domain == '2D' else trajectory[:, :3]
    min_coords = pos_data.min(axis=0) - args.padding
    max_coords = pos_data.max(axis=0) + args.padding
    grid_x = np.linspace(min_coords[0], max_coords[0], args.resolution)
    grid_y = np.linspace(min_coords[1], max_coords[1], args.resolution)
    grid_z = np.linspace(min_coords[2], max_coords[2], args.resolution) if sim_params.domain == '3D' else None

    if args.mode == 'static':
        run_static_mode(args, sim_params, trajectory, grid_x, grid_y, grid_z)
    elif args.mode == 'animation':
        run_animation_mode(args, sim_params, trajectory, grid_x, grid_y, grid_z)


if __name__ == '__main__':
    main()
