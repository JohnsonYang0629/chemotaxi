import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import matplotlib.patches as patches
from read_input import read_input
from read_input import read_vertex_file
from numba import njit, prange
import numba as nb
import os
import sys

# --- Custom Colormap Optimized for Black Background ---
# Bright Cyan (Sink) -> Black (0) -> Bright Orange (Source)
hc_colors = ["#00FFFF", "black", "#FF5500"]
bbr_cmap = LinearSegmentedColormap.from_list("cyan_black_orange", hc_colors)


# --- Calculation Functions ---
@njit(fastmath=True)
def green_function_2d(r_vec_sq, tau, pe, sigma=1.0):
    if tau <= 1e-12: return 0.0
    return (sigma * pe / (4 * np.pi * tau)) * np.exp(-r_vec_sq * pe / (4 * tau))


@njit(fastmath=True)
def green_function_3d(r_vec_sq, tau, pe, sigma=1.0):
    if tau <= 1e-12: return 0.0
    return (sigma * (pe / (4 * np.pi * tau)) ** 1.5) * np.exp(-r_vec_sq * pe / (4 * tau))


@njit(parallel=True, fastmath=True)
def _precompute_rotated_nodes_2d(n_history_steps, n_nodes, trajectory_history, nodes_body_frame):
    rotated_nodes_history = np.empty((n_history_steps, n_nodes, 2))
    for i in prange(n_history_steps):
        pos_center_hist = trajectory_history[i, :2]
        orient_hist = trajectory_history[i, 2:4]
        cos_t, sin_t = orient_hist[0], orient_hist[1]
        rot_matrix_hist = np.array([[cos_t, -sin_t], [sin_t, cos_t]])
        rotated_nodes_history[i, :, :] = pos_center_hist + (nodes_body_frame @ rot_matrix_hist.T)
    return rotated_nodes_history


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
def _precompute_rotated_nodes_3d(n_history_steps, n_nodes, trajectory_history, nodes_body_frame):
    rotated_nodes_history = np.empty((n_history_steps, n_nodes, 3))
    for i in prange(n_history_steps):
        pos_center_hist = trajectory_history[i, :3]
        orient_hist = trajectory_history[i, 3:7]
        rot_matrix_hist = quaternion_to_rotation_matrix(orient_hist)
        rotated_nodes_history[i, :, :] = pos_center_hist + (nodes_body_frame @ rot_matrix_hist.T)
    return rotated_nodes_history


@njit(parallel=True, fastmath=True)
def calculate_concentration_field_2d_unified(grid_points_x, grid_points_y, trajectories,
                                             nodes_body_frames_list, surface_sources_list,
                                             particle_type_flags, t1, dt, peclet_numbers):
    ny = len(grid_points_y)
    nx = len(grid_points_x)
    concentration_field = np.zeros((ny, nx))
    num_particles = trajectories.shape[0]

    for p in range(num_particles):
        trajectory = trajectories[p, :, :]
        pe = peclet_numbers[p]
        time_steps_to_consider = min(int(t1 / dt), trajectory.shape[0])

        if particle_type_flags[p] == 1:
            nodes_body_frame = nodes_body_frames_list[p]
            surface_sources = surface_sources_list[p]
            n_nodes = len(nodes_body_frame)

            rotated_nodes = _precompute_rotated_nodes_2d(
                time_steps_to_consider, n_nodes,
                trajectory[:time_steps_to_consider], nodes_body_frame
            )

            for i in prange(ny):
                y_val = grid_points_y[i]
                for j in range(nx):
                    x_val = grid_points_x[j]
                    particle_concentration = 0.0

                    for step in range(time_steps_to_consider):
                        tau = t1 - (step * dt)
                        if tau <= 1e-12:
                            continue

                        for s in range(n_nodes):
                            world_x = rotated_nodes[step, s, 0]
                            world_y = rotated_nodes[step, s, 1]
                            sigma = surface_sources[s]

                            r_vec_sq = (x_val - world_x) ** 2 + (y_val - world_y) ** 2
                            particle_concentration += green_function_2d(r_vec_sq, tau, pe, sigma) * dt

                    concentration_field[i, j] += particle_concentration

        else:
            for i in prange(ny):
                y_val = grid_points_y[i]
                for j in range(nx):
                    x_val = grid_points_x[j]
                    particle_concentration = 0.0

                    for step in range(time_steps_to_consider):
                        source_x = trajectory[step, 0]
                        source_y = trajectory[step, 1]
                        tau = t1 - (step * dt)
                        if tau <= 1e-12:
                            continue

                        r_vec_sq = (x_val - source_x) ** 2 + (y_val - source_y) ** 2
                        particle_concentration += green_function_2d(r_vec_sq, tau, pe, sigma=1.0) * dt

                    concentration_field[i, j] += particle_concentration

    return concentration_field


@njit(parallel=True)
def calculate_concentration_on_slices_3d_unified(grid_x, grid_y, grid_z, x_slice_idx, y_slice_idx, trajectories,
                                                 nodes_body_frames_list, surface_sources_list,
                                                 particle_type_flags, t1, dt, peclet_numbers):
    num_particles = trajectories.shape[0]
    yz_slice = np.zeros((len(grid_y), len(grid_z)))
    xz_slice = np.zeros((len(grid_x), len(grid_z)))

    x_val = grid_x[x_slice_idx]
    for j in prange(len(grid_y)):
        for k in prange(len(grid_z)):
            target_point = np.array([x_val, grid_y[j], grid_z[k]])
            total_concentration = 0.0
            for p in range(num_particles):
                trajectory = trajectories[p, :, :]
                pe = peclet_numbers[p]
                time_steps_to_consider = min(int(t1 / dt), trajectory.shape[0])
                particle_concentration = 0.0

                if particle_type_flags[p] == 1:
                    nodes_body_frame = nodes_body_frames_list[p]
                    surface_sources = surface_sources_list[p]
                    n_nodes = len(nodes_body_frame)
                    rotated_nodes = _precompute_rotated_nodes_3d(time_steps_to_consider, n_nodes,
                                                                 trajectory[:time_steps_to_consider], nodes_body_frame)
                    for step in range(time_steps_to_consider):
                        tau = t1 - (step * dt)
                        for s in range(n_nodes):
                            world_pos = rotated_nodes[step, s, :]
                            sigma = surface_sources[s]
                            r_vec_sq = np.sum((target_point - world_pos) ** 2)
                            particle_concentration += green_function_3d(r_vec_sq, tau, pe, sigma) * dt
                else:
                    for step in range(time_steps_to_consider):
                        source_pos = trajectory[step, :3]
                        tau = t1 - (step * dt)
                        r_vec_sq = np.sum((target_point - source_pos) ** 2)
                        particle_concentration += green_function_3d(r_vec_sq, tau, pe, sigma=1.0) * dt
                total_concentration += particle_concentration
            yz_slice[j, k] = total_concentration

    y_val = grid_y[y_slice_idx]
    for i in prange(len(grid_x)):
        for k in prange(len(grid_z)):
            target_point = np.array([grid_x[i], y_val, grid_z[k]])
            total_concentration = 0.0
            for p in range(num_particles):
                trajectory = trajectories[p, :, :]
                pe = peclet_numbers[p]
                time_steps_to_consider = min(int(t1 / dt), trajectory.shape[0])
                particle_concentration = 0.0

                if particle_type_flags[p] == 1:
                    nodes_body_frame = nodes_body_frames_list[p]
                    surface_sources = surface_sources_list[p]
                    n_nodes = len(nodes_body_frame)
                    rotated_nodes = _precompute_rotated_nodes_3d(time_steps_to_consider, n_nodes,
                                                                 trajectory[:time_steps_to_consider], nodes_body_frame)
                    for step in range(time_steps_to_consider):
                        tau = t1 - (step * dt)
                        for s in range(n_nodes):
                            world_pos = rotated_nodes[step, s, :]
                            sigma = surface_sources[s]
                            r_vec_sq = np.sum((target_point - world_pos) ** 2)
                            particle_concentration += green_function_3d(r_vec_sq, tau, pe, sigma) * dt
                else:
                    for step in range(time_steps_to_consider):
                        source_pos = trajectory[step, :3]
                        tau = t1 - (step * dt)
                        r_vec_sq = np.sum((target_point - source_pos) ** 2)
                        particle_concentration += green_function_3d(r_vec_sq, tau, pe, sigma=1.0) * dt
                total_concentration += particle_concentration
            xz_slice[i, k] = total_concentration

    return yz_slice.T, xz_slice.T


# --- Structure Loader Helper ---
def load_particle_structures(sim_params):
    nodes_body_frames_list = nb.typed.List()
    surface_sources_list = nb.typed.List()
    mobility_sources_list = nb.typed.List()
    particle_type_flags = []

    for i in range(sim_params.droplet_num):
        particle_type = sim_params.particle_types[i]
        structure_file = sim_params.structures[i][0] if i < len(sim_params.structures) else sim_params.structures[0][0]

        try:
            nodes = read_vertex_file.read_vertex_file(structure_file)
            nodes_body_frames_list.append(nodes)
        except Exception as e:
            print(f"FATAL: Could not load structure file for particle {i} from '{structure_file}'. Error: {e}",
                  file=sys.stderr)
            sys.exit(1)

        n_nodes = len(nodes)

        if hasattr(sim_params, 'mobility_distribution_files') and i < len(sim_params.mobility_distribution_files) and \
                sim_params.mobility_distribution_files[i] != 'None':
            mob_file = sim_params.mobility_distribution_files[i]
            try:
                mob = np.loadtxt(mob_file).flatten()
                mobility_sources_list.append(mob)
                print(f"Loaded Mobility data for particle {i} from '{mob_file}'.")
            except Exception as e:
                default_mob = sim_params.mobility_alphas[i]
                print(f"Warning: Could not load mobility file '{mob_file}': {e}. Using constant {default_mob}.")
                mobility_sources_list.append(np.full(n_nodes, default_mob, dtype=np.float64))
        else:
            default_mob = sim_params.mobility_alphas[i]
            mobility_sources_list.append(np.full(n_nodes, default_mob, dtype=np.float64))

        if particle_type == 'janus':
            particle_type_flags.append(1)
            chem_file = None
            if i < len(sim_params.chemical_distribution_files):
                chem_file = sim_params.chemical_distribution_files[i]
            elif len(sim_params.chemical_distribution_files) > 0:
                chem_file = sim_params.chemical_distribution_files[0]
                print(f"Info: Reusing chemical distribution from the first Janus particle for particle {i}.")

            if chem_file and chem_file != 'None':
                try:
                    sources = np.loadtxt(chem_file).flatten()
                    n_nodes = len(nodes_body_frames_list[-1])
                    sources = sources / n_nodes

                    surface_sources_list.append(sources)
                    print(f"Loaded Janus data for particle {i} from '{chem_file}'.")
                except Exception as e:
                    print(
                        f"FATAL: Could not load chemical distribution for Janus particle {i} from '{chem_file}'. Error: {e}",
                        file=sys.stderr)
                    sys.exit(1)
            else:
                print(f"FATAL: Particle {i} is 'janus' but no chemical_distribution_file was found.", file=sys.stderr)
                sys.exit(1)

        else:
            particle_type_flags.append(0)
            surface_sources_list.append(np.array([0.0]))
            print(f"Particle {i} is non-Janus.")

    return nodes_body_frames_list, surface_sources_list, mobility_sources_list, np.array(particle_type_flags,
                                                                                         dtype=np.int32)


# --- Particle Region Extraction Helper ---
def extract_signed_regions(angles_deg, vals):
    sort_idx = np.argsort(angles_deg)
    sorted_angles = angles_deg[sort_idx]
    sorted_vals = vals[sort_idx]

    max_val = np.max(np.abs(sorted_vals))
    threshold = 0.3 * max_val if max_val > 1e-5 else 1e-8

    states = np.zeros(len(sorted_vals), dtype=int)
    states[sorted_vals > threshold] = 1
    states[sorted_vals < -threshold] = -1

    regions = []
    if len(states) == 0 or max_val < 1e-5:
        return regions

    current_state = states[0]
    start_idx = 0
    n = len(sorted_angles)

    for i in range(1, n + 1):
        state = states[i % n] if i < n else states[0]

        if i == n or state != current_state:
            prev_idx = (start_idx - 1) % n
            diff_start = (sorted_angles[start_idx] - sorted_angles[prev_idx]) % 360
            t_start = sorted_angles[start_idx] - diff_start / 2.0

            end_idx = i - 1
            next_idx = i % n
            diff_end = (sorted_angles[next_idx] - sorted_angles[end_idx]) % 360
            t_end = sorted_angles[end_idx] + diff_end / 2.0

            while t_end < t_start:
                t_end += 360

            if current_state != 0:
                regions.append([current_state, t_start, t_end])

            current_state = state
            start_idx = i

    if len(regions) > 1 and regions[0][0] == regions[-1][0]:
        regions[0][1] = regions[-1][1] - 360
        regions.pop()

    return regions


# --- Plotting Functions ---
def plot_3d_sliced_view(yz_slice_data, xz_slice_data, grid_x, grid_y, grid_z, x_slice_idx, y_slice_idx, trajectories,
                        t1, dt, output_file, title, hide_markers=False):
    print("Generating 3D sliced view...")
    fig = plt.figure(figsize=(22, 7))
    ax1 = fig.add_subplot(1, 3, 1)
    ax2 = fig.add_subplot(1, 3, 2)
    ax3 = fig.add_subplot(1, 3, 3, projection='3d')
    fig.suptitle(title, fontsize=16)

    # Use TwoSlopeNorm for 3D slices as well
    c_min = min(yz_slice_data.min(),
                xz_slice_data.min()) if yz_slice_data.size > 0 and xz_slice_data.size > 0 else -1e-6
    c_max = max(yz_slice_data.max(), xz_slice_data.max()) if yz_slice_data.size > 0 and xz_slice_data.size > 0 else 1e-6

    # Safety bounds to ensure TwoSlopeNorm functions properly
    if c_min >= 0: c_min = -1e-6
    if c_max <= 0: c_max = 1e-6

    norm = TwoSlopeNorm(vcenter=0.0, vmin=c_min, vmax=c_max)

    ax1.set_facecolor('black')
    ax1.pcolormesh(grid_y, grid_z, yz_slice_data, cmap=bbr_cmap, norm=norm, shading='gouraud')
    ax1.set_title(f'YZ Cross-section at X = {grid_x[x_slice_idx]:.2f}')
    ax1.set_xlabel('Y')
    ax1.set_ylabel('Z')
    ax1.set_aspect('equal', 'box')

    ax2.set_facecolor('black')
    im = ax2.pcolormesh(grid_x, grid_z, xz_slice_data, cmap=bbr_cmap, norm=norm, shading='gouraud')
    ax2.set_title(f'XZ Cross-section at Y = {grid_y[y_slice_idx]:.2f}')
    ax2.set_xlabel('X')
    ax2.set_ylabel('Z')
    ax2.set_aspect('equal', 'box')
    fig.colorbar(im, ax=ax2, orientation='vertical', fraction=0.046, pad=0.04, label='Concentration C')

    all_traj_points = np.vstack(
        [trajectories[i, :min(int(t1 / dt), trajectories.shape[1]), :3] for i in range(trajectories.shape[0])])
    center_point = np.mean(all_traj_points, axis=0)

    for p_idx in range(trajectories.shape[0]):
        trajectory = trajectories[p_idx, :, :]
        time_steps_to_plot = min(int(t1 / dt), trajectory.shape[0])
        traj_to_plot = trajectory[:time_steps_to_plot, :3]
        ax1.plot(traj_to_plot[:, 1], traj_to_plot[:, 2], 'w-', lw=1.5, alpha=0.7)
        ax2.plot(traj_to_plot[:, 0], traj_to_plot[:, 2], 'w-', lw=1.5, alpha=0.7)

        if not hide_markers:
            ax1.scatter(traj_to_plot[-1, 1], traj_to_plot[-1, 2], c='red', s=40, zorder=5)
            ax2.scatter(traj_to_plot[-1, 0], traj_to_plot[-1, 2], c='red', s=40, zorder=5)

        points = traj_to_plot.reshape(-1, 1, 3)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        t_values = np.linspace(0, t1, len(traj_to_plot))
        lc = Line3DCollection(segments, cmap=plt.get_cmap('plasma'),
                              norm=Normalize(vmin=t_values.min(), vmax=t_values.max()))
        lc.set_array(t_values)
        lc.set_linewidth(2)
        ax3.add_collection(lc)

        if not hide_markers:
            start_label = 'Start' if p_idx == 0 else None
            end_label = f'End (t={t1:.1f})' if p_idx == 0 else None
            ax3.scatter(traj_to_plot[0, 0], traj_to_plot[0, 1], traj_to_plot[0, 2], c='lime', s=50, label=start_label,
                        depthshade=False)
            ax3.scatter(traj_to_plot[-1, 0], traj_to_plot[-1, 1], traj_to_plot[-1, 2], c='red', s=50, label=end_label,
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

    handles, labels = ax3.get_legend_handles_labels()
    if handles and not hide_markers:
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

def run_static_mode(args, sim_params, trajectories, grid_x, grid_y, grid_z):
    print(f"--- Running in STATIC mode for time t = {args.t1} ---")
    output_filename = args.output_file if (args.output_file.lower().endswith('.png') or
                                           args.output_file.lower().endswith('.gif') or
                                           args.output_file.lower().endswith('.mp4')) else args.output_file + '.png'

    nodes_body_frames_list, surface_sources_list, mobility_sources_list, particle_type_flags = load_particle_structures(
        sim_params)

    if sim_params.domain == '2D':
        print("Calculating 2D concentration field...")
        concentration_field = calculate_concentration_field_2d_unified(
            grid_x, grid_y, trajectories, nodes_body_frames_list, surface_sources_list,
            particle_type_flags, args.t1, sim_params.dt, sim_params.peclet_numbers
        )
        print("Calculation complete.")
        title = f'Chemical Concentration at t={args.t1:.2f}'

        fig, ax = plt.subplots(figsize=(10, 8))
        ax.set_facecolor('black')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_title(title)
        ax.set_aspect('equal', 'box')

        # Advanced normalization for asymmetric source/sink ranges
        c_min = concentration_field.min()
        c_max = concentration_field.max()
        if c_min >= 0: c_min = -1e-6
        if c_max <= 0: c_max = 1e-6

        norm = TwoSlopeNorm(vcenter=0.0, vmin=c_min, vmax=c_max)

        im = ax.pcolormesh(grid_x, grid_y, concentration_field, shading='gouraud', cmap=bbr_cmap, norm=norm, zorder=1)
        ax.grid(False)
        fig.colorbar(im, ax=ax, label='Chemical Concentration C')

        # --- GRADIENT PLOTTING BLOCK ---
        if args.plot_gradient:
            print("Calculating gradient field...")
            dy = grid_y[1] - grid_y[0]
            dx = grid_x[1] - grid_x[0]
            grad_y, grad_x = np.gradient(concentration_field, dy, dx)

            # Calculate gradient magnitude and log magnitude (with zero-protection)
            grad_mag = np.sqrt(grad_x ** 2 + grad_y ** 2)
            safe_mag = np.where(grad_mag > 1e-12, grad_mag, 1e-12)

            # Normalize vectors to ensure uniform quiver arrow lengths
            U = grad_x / safe_mag
            V = grad_y / safe_mag
            ln_grad_mag = np.log(safe_mag)

            grid_X, grid_Y = np.meshgrid(grid_x, grid_y)
            skip = args.quiver_density

            # Using negative gradients (-U, -V) for physical flow/diffusion direction.
            # Blues_r gives bright cyan/white to high magnitudes, standing out on dark backgrounds.
            q = ax.quiver(grid_X[::skip, ::skip], grid_Y[::skip, ::skip],
                          -U[::skip, ::skip], -V[::skip, ::skip],
                          ln_grad_mag[::skip, ::skip],
                          cmap='Blues_r', pivot='mid', scale=40,
                          width=0.003, headwidth=4, headlength=4, zorder=2)

            cbar_grad = fig.colorbar(q, ax=ax, orientation='horizontal', fraction=0.046, pad=0.1)
            cbar_grad.set_label(r'$\ln|-\nabla C|$', fontsize=12)
        # -------------------------------

        for i in range(trajectories.shape[0]):
            trajectory = trajectories[i, :, :]
            time_steps_to_plot = min(int(args.t1 / sim_params.dt), trajectory.shape[0])
            traj_to_plot = trajectory[:time_steps_to_plot]
            ax.plot(traj_to_plot[:, 0], traj_to_plot[:, 1], 'w-', lw=1.5, alpha=0.8, zorder=2)

            if not args.hide_markers:
                ax.plot(traj_to_plot[0, 0], traj_to_plot[0, 1], 'go', markersize=8, label='Start', zorder=3)
                ax.plot(traj_to_plot[-1, 0], traj_to_plot[-1, 1], 'ro', markersize=8,
                        label=f'Current (t={args.t1:.2f})',
                        zorder=3)

            if args.visualize_particle != 'none':
                nodes = nodes_body_frames_list[i]
                radius = np.max(np.sqrt(np.sum(nodes[:, :2] ** 2, axis=1)))
                current_x, current_y = traj_to_plot[-1, 0], traj_to_plot[-1, 1]

                circle = patches.Circle((current_x, current_y), radius, facecolor='white', edgecolor='none', alpha=0.5,
                                        zorder=4)
                ax.add_patch(circle)

                if args.visualize_particle in ['chem', 'mobility', 'all']:
                    cos_t, sin_t = traj_to_plot[-1, 2], traj_to_plot[-1, 3]
                    theta_rot_deg = np.degrees(np.arctan2(sin_t, cos_t))

                    angles_rad = np.arctan2(nodes[:, 1], nodes[:, 0])
                    angles_deg = np.degrees(angles_rad)

                    if args.visualize_particle in ['chem', 'all']:
                        regions = extract_signed_regions(angles_deg, surface_sources_list[i])
                        for state, t_start, t_end in regions:
                            color = 'lightcoral' if state == 1 else 'cornflowerblue'
                            w = patches.Wedge((current_x, current_y), radius,
                                              t_start + theta_rot_deg,
                                              t_end + theta_rot_deg,
                                              facecolor=color, edgecolor='none', alpha=0.5, zorder=5)
                            ax.add_patch(w)

                    if args.visualize_particle in ['mobility', 'all']:
                        regions = extract_signed_regions(angles_deg, mobility_sources_list[i])
                        for state, t_start, t_end in regions:
                            color = 'cadetblue' if state == 1 else 'palevioletred'
                            w = patches.Wedge((current_x, current_y), radius,
                                              t_start + theta_rot_deg,
                                              t_end + theta_rot_deg,
                                              width=radius * 0.15,
                                              facecolor=color, edgecolor='none', alpha=0.8, zorder=6)
                            ax.add_patch(w)

        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))

        if args.range:
            ax.set_xlim(args.range[0], args.range[1])
            ax.set_ylim(args.range[2], args.range[3])
        elif args.center:
            ax.set_xlim(grid_x.min(), grid_x.max())
            ax.set_ylim(grid_y.min(), grid_y.max())

        plt.savefig(output_filename, dpi=300)
        print(f"Static 2D plot saved to {output_filename}")

    elif sim_params.domain == '3D':
        print("Preparing data for 3D concentration calculation...")

        if args.slice_center == 'auto':
            all_traj_points = np.vstack(
                [trajectories[i, :min(int(args.t1 / sim_params.dt), trajectories.shape[1]), :3] for i in
                 range(trajectories.shape[0])])
            center_point = np.mean(all_traj_points, axis=0)
            print(f"Slice center method: 'auto'. Calculated center: {center_point}")
        elif args.slice_center == 'end':
            end_step = min(int(args.t1 / sim_params.dt), trajectories.shape[1] - 1)
            center_point = trajectories[0, end_step, :3]
            print(f"Slice center method: 'end'. Using center: {center_point}")
        else:
            try:
                center_point = np.array(list(map(float, args.slice_center.split(','))))
                if center_point.shape != (3,):
                    raise ValueError
                print(f"Slice center method: 'manual'. Using specified center: {center_point}")
            except (ValueError, AttributeError):
                print(
                    f"Error: Invalid format for --slice-center. Use 'auto', 'end', or 'x,y,z'. You provided: '{args.slice_center}'",
                    file=sys.stderr)
                sys.exit(1)

        x_slice_idx = np.argmin(np.abs(grid_x - center_point[0]))
        y_slice_idx = np.argmin(np.abs(grid_y - center_point[1]))

        print("Calculating concentration on 2D slices (this may take a moment)...")
        yz_slice, xz_slice = calculate_concentration_on_slices_3d_unified(
            grid_x, grid_y, grid_z, x_slice_idx, y_slice_idx, trajectories,
            nodes_body_frames_list, surface_sources_list, particle_type_flags,
            args.t1, sim_params.dt, sim_params.peclet_numbers
        )

        print("Calculation complete.")
        title = f'3D Sliced View at t={args.t1:.2f}'
        plot_3d_sliced_view(yz_slice, xz_slice, grid_x, grid_y, grid_z, x_slice_idx, y_slice_idx, trajectories, args.t1,
                            sim_params.dt,
                            output_filename, title, args.hide_markers)


def run_animation_mode(args, sim_params, trajectories, grid_x, grid_y, grid_z):
    if sim_params.domain == '3D':
        print("Animation mode for 3D is not yet implemented. Please use static mode.", file=sys.stderr)
        return

    total_time = trajectories.shape[1] * sim_params.dt
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

    nodes_body_frames_list, surface_sources_list, mobility_sources_list, particle_type_flags = load_particle_structures(
        sim_params)

    print("Pre-calculating color limits for stable animation...")
    final_field = calculate_concentration_field_2d_unified(
        grid_x, grid_y, trajectories, nodes_body_frames_list, surface_sources_list,
        particle_type_flags, total_time, sim_params.dt, sim_params.peclet_numbers
    )

    # Asymmetric mapping for animations
    vmin, vmax = final_field.min(), final_field.max()
    if vmin >= 0: vmin = -1e-6
    if vmax <= 0: vmax = 1e-6
    norm = TwoSlopeNorm(vcenter=0.0, vmin=vmin, vmax=vmax)

    print(f"Global color range set via TwoSlopeNorm: [{vmin:.4f}, {vmax:.4f}]")

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.set_facecolor('black')
    title_obj = ax.set_title('')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_aspect('equal', 'box')

    im = ax.pcolormesh(grid_x, grid_y, np.zeros_like(final_field),
                       shading='gouraud',
                       cmap=bbr_cmap,
                       norm=norm,
                       zorder=1)

    lines = [ax.plot([], [], 'w-', lw=1.5, alpha=0.8, zorder=2)[0] for _ in range(trajectories.shape[0])]
    start_dots = [ax.plot(trajectories[i, 0, 0], trajectories[i, 0, 1], 'go', markersize=8, zorder=3)[0] for i in
                  range(trajectories.shape[0])]
    current_dots = [ax.plot([], [], 'ro', markersize=8, zorder=3)[0] for _ in range(trajectories.shape[0])]

    particle_circles = []
    particle_chem_dots = []
    particle_mob_dots = []

    if args.visualize_particle != 'none':
        for i in range(trajectories.shape[0]):
            nodes = nodes_body_frames_list[i]
            radius = np.max(np.sqrt(np.sum(nodes[:, :2] ** 2, axis=1)))

            start_x, start_y = trajectories[i, 0, 0], trajectories[i, 0, 1]
            circle = patches.Circle((start_x, start_y), radius, facecolor='white', edgecolor='none', alpha=0.5,
                                    zorder=4)
            ax.add_patch(circle)
            particle_circles.append(circle)

            cos_t, sin_t = trajectories[i, 0, 2], trajectories[i, 0, 3]
            theta_rot_deg = np.degrees(np.arctan2(sin_t, cos_t))

            angles_rad = np.arctan2(nodes[:, 1], nodes[:, 0])
            angles_deg = np.degrees(angles_rad)

            if args.visualize_particle in ['chem', 'all']:
                regions = extract_signed_regions(angles_deg, surface_sources_list[i])
                chem_wedges = []
                for state, t_start, t_end in regions:
                    color = 'lightcoral' if state == 1 else 'cornflowerblue'
                    w = patches.Wedge((start_x, start_y), radius,
                                      t_start + theta_rot_deg, t_end + theta_rot_deg,
                                      facecolor=color, edgecolor='none', alpha=0.5, zorder=5)
                    ax.add_patch(w)
                    chem_wedges.append((w, t_start, t_end))
                particle_chem_dots.append(chem_wedges)
            else:
                particle_chem_dots.append(None)

            if args.visualize_particle in ['mobility', 'all']:
                regions = extract_signed_regions(angles_deg, mobility_sources_list[i])
                mob_wedges = []
                for state, t_start, t_end in regions:
                    color = 'cadetblue' if state == 1 else 'palevioletred'
                    w = patches.Wedge((start_x, start_y), radius,
                                      t_start + theta_rot_deg, t_end + theta_rot_deg,
                                      width=radius * 0.15,
                                      facecolor=color, edgecolor='none', alpha=0.8, zorder=6)
                    ax.add_patch(w)
                    mob_wedges.append((w, t_start, t_end))
                particle_mob_dots.append(mob_wedges)
            else:
                particle_mob_dots.append(None)

    ax.legend([start_dots[0], current_dots[0]], ['Start', 'Current Position'], loc='upper right', facecolor='white',
              framealpha=0.8)
    ax.grid(False)

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label('Chemical Concentration C')

    def update(frame_time):
        concentration_field = calculate_concentration_field_2d_unified(
            grid_x, grid_y, trajectories, nodes_body_frames_list, surface_sources_list,
            particle_type_flags, frame_time, sim_params.dt, sim_params.peclet_numbers
        )
        im.set_array(concentration_field.ravel())
        render_objs = [im] + lines + current_dots + [title_obj]

        for i in range(trajectories.shape[0]):
            trajectory = trajectories[i, :, :]
            current_step = min(int(frame_time / sim_params.dt), trajectory.shape[0] - 1)
            traj_to_plot = trajectory[:current_step + 1]
            lines[i].set_data(traj_to_plot[:, 0], traj_to_plot[:, 1])
            current_dots[i].set_data(traj_to_plot[-1, 0], traj_to_plot[-1, 1])

            if args.visualize_particle != 'none':
                current_x, current_y = traj_to_plot[-1, 0], traj_to_plot[-1, 1]
                particle_circles[i].center = (current_x, current_y)
                render_objs.append(particle_circles[i])

                if args.visualize_particle in ['chem', 'mobility', 'all']:
                    cos_t, sin_t = traj_to_plot[-1, 2], traj_to_plot[-1, 3]
                    theta_rot_deg = np.degrees(np.arctan2(sin_t, cos_t))

                    if particle_chem_dots[i] is not None:
                        for w, t_start, t_end in particle_chem_dots[i]:
                            w.set_center((current_x, current_y))
                            w.set_theta1(t_start + theta_rot_deg)
                            w.set_theta2(t_end + theta_rot_deg)
                            render_objs.append(w)

                    if particle_mob_dots[i] is not None:
                        for w, t_start, t_end in particle_mob_dots[i]:
                            w.set_center((current_x, current_y))
                            w.set_theta1(t_start + theta_rot_deg)
                            w.set_theta2(t_end + theta_rot_deg)
                            render_objs.append(w)

        title_obj.set_text(f't={frame_time:.2f} s')
        print(f"Rendering frame for time: {frame_time:.2f} / {total_time:.2f}")

        return render_objs

    print("Starting animation rendering... This may take a while.")
    try:
        ani = animation.FuncAnimation(fig, update, frames=animation_times, blit=False, repeat=False)
        ani.save(output_filename, writer=writer_instance)
        print(f"Animation successfully saved to {output_filename}")
    except FileNotFoundError:
        print(f"\n--- ERROR: '{args.writer}' writer not found ---", file=sys.stderr)
        if args.writer == 'ffmpeg':
            print("To save MP4 videos, please install ffmpeg.", file=sys.stderr)
        else:
            print("To save GIFs, please install the Pillow library: pip install Pillow", file=sys.stderr)
        sys.exit(1)


def read_trajectories(filepath, num_particles):
    trajectories = [[] for _ in range(num_particles)]
    with open(filepath, 'r') as f:
        lines = f.readlines()

    lines_per_frame = num_particles + 1
    num_frames = len(lines) // lines_per_frame

    for frame in range(num_frames):
        start_index = frame * lines_per_frame

        for p in range(num_particles):
            line_index = start_index + 1 + p
            line = lines[line_index]
            parts = list(map(float, line.split()))
            trajectories[p].append(parts)

    return np.array(trajectories)


def main():
    parser = argparse.ArgumentParser(
        description="Visualize chemical concentration from a chemotaxi simulation statically or as an animation.")

    parser.add_argument('--file-prefix', dest='file_prefix', type=str, default=None,
                        help='Path and prefix for input files. If used, --input-file and --trajectory-file are ignored.')

    parser.add_argument('--input-file', dest='input_file', type=str,
                        help='(Optional) Full path to simulation input file. Use if prefix is not provided.')
    parser.add_argument('--trajectory-file', dest='trajectory_file', type=str,
                        help='(Optional) Full path to trajectory file. Use if prefix is not provided.')

    parser.add_argument('--output-file', dest='output_file', type=str, default=None,
                        help='Base name of output file (no extension). Defaults to prefix or "concentration_field".')
    parser.add_argument('--mode', choices=['static', 'animation'], default='static',
                        help="Choose visualization mode. Default: static.")
    parser.add_argument('--writer', choices=['ffmpeg', 'gif'], default='gif',
                        help="[ANIMATION mode] Writer to use for saving animation. 'ffmpeg' for MP4, 'gif' for GIF. Default: gif.")
    parser.add_argument('--time', dest='t1', type=float, help='[STATIC mode] Specific time t1 to visualize.')
    parser.add_argument('--slice-center', dest='slice_center', type=str, default='auto',
                        help="[STATIC 3D mode] Method to center slices: 'auto', 'end', or 'x,y,z'.")
    parser.add_argument('--frame-interval', type=float, default=5.0,
                        help='[ANIMATION mode] Simulation time between frames. Default: 5.0.')
    parser.add_argument('--resolution', type=int, default=100, help='Grid resolution for visualization.')
    parser.add_argument('--padding', type=float, default=5.0, help='Padding around the trajectory for the grid.')

    parser.add_argument('--visualize-particle', nargs='?', const='solid', default='none',
                        choices=['none', 'solid', 'chem', 'mobility', 'all'],
                        help='[2D mode] Visualize particle: none, solid, chem, mobility, or all (double rings).')

    parser.add_argument('--hide-markers', action='store_true',
                        help='[STATIC mode] Hide the start/current position markers and legend.')

    parser.add_argument('--center', type=float, nargs=2, metavar=('X0', 'Y0'),
                        help='[2D STATIC mode] Center the visualization at (X0, Y0). Keeps the original width/height.')
    parser.add_argument('--range', type=float, nargs=4, metavar=('XMIN', 'XMAX', 'YMIN', 'YMAX'),
                        help='[2D STATIC mode] Explicitly set the visualization grid range. Overrides --center.')

    parser.add_argument('--plot-gradient', action='store_true',
                        help='[2D STATIC mode] Overlay the gradient field using arrows.')
    parser.add_argument('--quiver-density', type=int, default=4,
                        help='[2D STATIC mode] Density of gradient arrows (plots every Nth point). Default: 4.')

    args = parser.parse_args()

    if args.file_prefix:
        input_file = args.file_prefix + '.inputfile'
        trajectory_file = args.file_prefix + '.config'
        if args.output_file is None:
            args.output_file = args.file_prefix + '.concentration_field'

        if not os.path.exists(input_file) or not os.path.exists(trajectory_file):
            print(f"Error: Using prefix '{args.file_prefix}', but could not find one or both files:", file=sys.stderr)
            print(f"  - Attempted input file: {input_file}", file=sys.stderr)
            print(f"  - Attempted trajectory file: {trajectory_file}", file=sys.stderr)
            print("Please check the prefix or provide full paths using --input-file and --trajectory-file.",
                  file=sys.stderr)
            sys.exit(1)

    elif args.input_file and args.trajectory_file:
        input_file = args.input_file
        trajectory_file = args.trajectory_file
        if args.output_file is None:
            args.output_file = 'concentration_field'
    else:
        parser.error("You must provide either --file-prefix or both --input-file and --trajectory-file.")

    if args.mode == 'static' and args.t1 is None:
        parser.error("--time is required for --mode='static'")

    print("Reading input files...")
    sim_params = read_input.ReadInput(input_file)
    trajectories = read_trajectories(trajectory_file, sim_params.droplet_num)

    print(f"Successfully loaded {trajectories.shape[0]} particle trajectories.")
    print(
        f"Simulation Parameters: Domain={sim_params.domain}, Peclet Numbers={sim_params.peclet_numbers}, dt={sim_params.dt}")

    all_pos_data = np.vstack([trajectories[i, :, :2] if sim_params.domain == '2D' else trajectories[i, :, :3] for i in
                              range(trajectories.shape[0])])
    min_coords = all_pos_data.min(axis=0) - args.padding
    max_coords = all_pos_data.max(axis=0) + args.padding

    if sim_params.domain == '2D':
        if args.range:
            min_coords[0], max_coords[0], min_coords[1], max_coords[1] = args.range
        elif args.center:
            cx, cy = args.center
            width = max_coords[0] - min_coords[0]
            height = max_coords[1] - min_coords[1]
            min_coords[0] = cx - width / 2.0
            max_coords[0] = cx + width / 2.0
            min_coords[1] = cy - height / 2.0
            max_coords[1] = cy + height / 2.0

    grid_x = np.linspace(min_coords[0], max_coords[0], args.resolution)
    grid_y = np.linspace(min_coords[1], max_coords[1], args.resolution)
    grid_z = np.linspace(min_coords[2], max_coords[2], args.resolution) if sim_params.domain == '3D' else None

    if args.mode == 'static':
        run_static_mode(args, sim_params, trajectories, grid_x, grid_y, grid_z)
    elif args.mode == 'animation':
        run_animation_mode(args, sim_params, trajectories, grid_x, grid_y, grid_z)


if __name__ == '__main__':
    main()
