import numpy as np
import scipy.special as sc
import sys
from scipy.special import kn, k0, k1, erfc
import math
from scipy import integrate
from numba import njit, prange
from numba import jit
import time
try:
    import multiprocessing as mp
except ImportError:
    print('Parallel module missing')


@jit(nopython=True, fastmath=True)
def calc_gradient_history_part_2d(target_position, history_path_location, peclet_number, step, dt):
    # Loop accelerated with numba
    # Recommended method
    # For 2D cases

    gradient_segment = np.zeros((step, 2))
    for s in range(step):
        position_difference = target_position - history_path_location[s, :]
        position_difference_norm_square = np.power(np.linalg.norm(position_difference), 2)
        gradient_segment[s, :] = 2 * np.power(peclet_number/4, 2)/(np.pi * np.power((step - s) * dt, 2)) \
                                   * np.exp(-peclet_number/(4 * (step - s) * dt) * position_difference_norm_square)\
                                   * position_difference
    # Use trapezoidal rule to sum up segments
    gradient_history_part = 1/2 * dt * np.sum(gradient_segment, axis=0)
    return gradient_history_part


@jit(nopython=True, fastmath=True)
def calc_gradient_history_part_2d_optimized(target_position, history_path_location, peclet_number, step, dt):
    gradient_segment = np.zeros((step, 2))
    for s in range(step):
        position_difference = target_position - history_path_location[s, :]

        position_difference_norm_square = position_difference[0] ** 2 + position_difference[1] ** 2

        tau = (step - s) * dt

        if tau <= 1e-12: continue

        coeff = 2 * (peclet_number / 4) ** 2 / (np.pi * tau ** 2)
        exp_term = np.exp(-peclet_number / (4 * tau) * position_difference_norm_square)
        gradient_segment[s, :] = coeff * exp_term * position_difference

    gradient_history_part = 0.5 * dt * np.sum(gradient_segment, axis=0)
    return gradient_history_part


@jit(nopython=True, fastmath=True)
def calc_gradient_local_part_2d_optimized(target_position, history_path_location, peclet_number, step, dt):
    position_difference = target_position - history_path_location[step - 1, :]

    position_difference_norm_square = position_difference[0] ** 2 + position_difference[1] ** 2

    # Small position differences lead to large local part
    if position_difference_norm_square < 1e-1:
        return np.zeros(2)

    rho = peclet_number * position_difference_norm_square / (4 * dt)

    common_factor = 2 * peclet_number / (4 * np.pi * position_difference_norm_square ** 2)

    gradient_local_part = common_factor * np.exp(-rho) * (rho ** 2 + 2 * rho + 2) * position_difference

    return gradient_local_part


@jit(nopython=True, fastmath=True)
def calc_gradient_history_part_3d(target_position, history_path_location, peclet_number, step, dt):
    # Loop accelerated with numba
    # Recommended method
    # For 3D cases

    gradient_segment = np.zeros((step, 3))
    for s in range(step):
        position_difference = target_position - history_path_location[s, :]
        position_difference_norm_square = np.power(np.linalg.norm(position_difference), 2)
        gradient_segment[s, :] = 2 * np.power(np.pi/4, 3/2) \
                                   * np.power(peclet_number/4, 5/2) / (np.power((step - s) * dt, 5/2)) \
                                   * np.exp(-peclet_number/(4 * (step - s) * dt) * position_difference_norm_square)\
                                   * position_difference
    # Use trapezoidal rule to sum up segments
    gradient_history_part = 1/2 * dt * np.sum(gradient_segment, axis=0)
    return gradient_history_part


@jit(nopython=True, fastmath=True)
def calc_gradient_history_part_np_acc_2d(target_position, history_path_location, peclet_number, step, dt):
    # Accelerated with numpy matrix operation
    # Computational time larger than numba accelerated loop
    # NOT recommended

    gradient_segment = np.zeros((step, 2))
    target_position = target_position * np.ones((step, 2))
    history_path_location_current = history_path_location[0:step, :]
    position_difference = target_position - history_path_location_current
    position_difference_norm_square = np.sum(np.abs(position_difference)**2, axis=-1)
    time_step = dt * np.linspace(step, 1, step)
    gradient_segment_coeff = 2 * np.power(peclet_number / 4, 2) / (np.pi * np.power(time_step, 2)) \
                               * np.exp(-peclet_number / (4 * time_step) * position_difference_norm_square)
    gradient_segment = np.dot(np.diag(gradient_segment_coeff), position_difference)
    # Use trapezoidal rule to sum up segments
    gradient_history_part = 1/2 * dt * np.sum(gradient_segment, axis=0)
    return gradient_history_part


def calc_gradient_history_part_parallel_2d(target_position, history_path_location, peclet_number, step, dt):
    # Accelerated with Python built-in parallel computing
    # Parallel Pools: (N-1) steps
    # Potentially useful for EXTRA-long simulation

    step_list = np.arange(step)
    num_cores = int(mp.cpu_count())
    process_pool = mp.Pool(num_cores - 2)
    start_time = time.time()
    gradient_segment_results = [process_pool.apply_async(calc_gradient_history_segment_2d,
                                args=(target_position, history_path_location, peclet_number, step, current_step, dt))
                                for current_step in step_list]
    process_pool.close()
    gradient_segment = [segment.get() for segment in gradient_segment_results]
    elapsed_time = time.time() - start_time
    print('parallel process time = ', elapsed_time)
    # Use trapezoidal rule to sum up segments
    gradient_history_part = 1/2 * dt * np.sum(gradient_segment)
    return gradient_history_part


def calc_gradient_2d_parallel(seq, target_points_abs_loc, location_history, structure_ref_config, peclet_number, step, dt):
    # Complimentary for function: calc_surface_gradient_circle

    gradient_history_part = calc_gradient_history_part_2d(target_points_abs_loc[seq, :],
                                                          location_history, peclet_number, step, dt)
    gradient_local_part = calc_gradient_local_part_2d(target_points_abs_loc[seq, :],
                                                      location_history, peclet_number, step, dt)
    surface_gradient = calc_tangential_gradient_part(structure_ref_config, seq,
                                                     gradient_history_part + gradient_local_part)
    return surface_gradient


def calc_gradient_3d_parallel(seq, target_points_abs_loc, location_history, structure_ref_config, peclet_number, step, dt):
    # Complimentary for function: calc_surface_gradient_circle

    gradient_history_part = calc_gradient_history_part_3d(target_points_abs_loc[seq, :],
                                                          location_history, peclet_number, step, dt)
    #gradient_local_part = calc_gradient_local_part_2d(target_points_abs_loc[seq, :],
    #                                                  location_history, peclet_number, step, dt)
    surface_gradient = calc_tangential_gradient_part(structure_ref_config, seq,
                                                     gradient_history_part)  # + gradient_local_part)
    return surface_gradient


def calc_gradient_history_segment_2d(target_position, history_path_location, peclet_number, step, current_step, dt):
    # Complimentary part for "calc_gradient_history_part_parallel_2d"
    # Calculate the chemical gradient for ith step in history part [0, (N-1)dt]
    position_difference = target_position - history_path_location[current_step, :]
    position_difference_norm_square = np.power(np.linalg.norm(position_difference), 2)
    gradient_segment = 2 * np.power(peclet_number / 4, 2) / (np.pi * np.power((step - current_step) * dt, 2)) \
                         * np.exp(-peclet_number / (4 * (step - current_step) * dt) * position_difference_norm_square) \
                         * position_difference
    return gradient_segment


def calc_gradient_local_part_2d(target_position, history_path_location, peclet_number, step, dt):
    # Calculate the chemical gradient for Nth step in local part [(N-1)dt, Ndt]
    # Complimentary for function: calc_surface_gradient_circle and calc_gradient_2d_parallel
    position_difference = target_position - history_path_location[step - 1, :]
    position_difference_norm_square = np.power(np.linalg.norm(position_difference), 2)
    rho = peclet_number * position_difference_norm_square/(4 * dt)
    gradient_local_part = 2 * peclet_number/(4 * np.pi * np.power(position_difference_norm_square, 2)) \
                            * np.exp(-rho) * (np.power(rho, 2) + 2 * rho + 2) \
                            * position_difference
    return gradient_local_part


def calc_tangential_gradient_part(structure_ref_config, seq, gradient):
    # Find the tangential part of the chemical gradient on a circle
    # 2D cases
    pos = structure_ref_config[seq, :]
    pos_norm = pos/np.linalg.norm(pos)
    # Given that discretized points are on an unit circle
    pos_tangent_vec = np.array([pos_norm[1], -pos_norm[0]])
    tangential_gradient_part = np.dot(gradient, pos_tangent_vec) * pos_tangent_vec
    return tangential_gradient_part


@njit(fastmath=True)
def calc_tangential_gradient_part_numba(structure_ref_config, seq, chemical_gradient_sum):
    """
    Calculates the tangential component of a gradient at a specific point on the structure (Numba-optimized).
    """
    # 1. Get the reference position vector for the current point
    pos_ref = structure_ref_config[seq]

    # 2. Normalize it to get the unit radial vector
    pos_ref_norm = np.sqrt(pos_ref[0] ** 2 + pos_ref[1] ** 2)
    if pos_ref_norm < 1e-12:
        return np.zeros(2)

    radial_vec_x = pos_ref[0] / pos_ref_norm
    radial_vec_y = pos_ref[1] / pos_ref_norm

    # 3. Calculate the unit tangential vector (by rotating the radial vector 90 degrees)
    tangent_vec_x = -radial_vec_y
    tangent_vec_y = radial_vec_x

    # 4. Project the chemical gradient onto the tangential vector
    # Projection = (gradient ⋅ tangent_vector) * tangent_vector
    dot_product = chemical_gradient_sum[0] * tangent_vec_x + chemical_gradient_sum[1] * tangent_vec_y

    tangential_gradient_x = dot_product * tangent_vec_x
    tangential_gradient_y = dot_product * tangent_vec_y

    return np.array([tangential_gradient_x, tangential_gradient_y])


def calc_tangential_gradient_part_3d(structure_ref_config, seq, gradient):
    # Find the tangential part of the chemical gradient on a sphere
    # 3D cases
    pos = structure_ref_config[seq, :]
    pos_norm = pos/np.linalg.norm(pos)
    # Given that discretized points are on a unit sphere
    pos_sph = cart2sph_vector_3d(pos_norm)
    pos_tangent_vec = sph2cart_field_3d(pos_sph, np.array([[0, 1, 0]]))
    tangential_gradient_part = np.dot(gradient, pos_tangent_vec) * pos_tangent_vec
    return tangential_gradient_part


def calc_surface_gradient_circle(body, peclet_number, structure_ref_config, dt, *args, **kwargs):
    chem_gradient = np.zeros([2])
    step = kwargs.get('step')
    location = body.location
    location_history = body.location_history
    target_points_abs_loc = structure_ref_config + location
    acceleration = kwargs.get('acceleration')
    num_core = kwargs.get('core')

    # First guess
    if step == 0 and step >= 0:
        chem_gradient = np.zeros([2])
    # After initialized: second step only have
    elif step == 1 and step >= 0:
        surface_gradient_history_sum = np.zeros([2])
        for seq, loc in enumerate(target_points_abs_loc):
            gradient_history_part = calc_gradient_history_part_2d(loc, location_history, peclet_number, step, dt)
            surface_gradient_history = calc_tangential_gradient_part(structure_ref_config, seq, gradient_history_part)
            surface_gradient_history_sum += surface_gradient_history

        chem_gradient = 2 * np.pi/np.shape(structure_ref_config)[0] * surface_gradient_history_sum

    elif step >= 2 and acceleration == "numba":
        # Numba accelerated
        # Averaged calculation time/step < 0.1s for step<=10000
        surface_gradient_sum = np.zeros([2])
        for seq, loc in enumerate(target_points_abs_loc):
            gradient_history_part = calc_gradient_history_part_2d(loc, location_history, peclet_number, step, dt)
            gradient_local_part = calc_gradient_local_part_2d(loc, location_history, peclet_number, step, dt)
            chemical_gradient_sum = gradient_history_part + gradient_local_part
            surface_gradient = calc_tangential_gradient_part(structure_ref_config, seq, chemical_gradient_sum)
            surface_gradient_sum += surface_gradient

        chem_gradient = 2 * np.pi/np.shape(structure_ref_config)[0] * surface_gradient_sum

    elif step >= 2 and acceleration == "numba_optimized":
        # Numba accelerated
        # Averaged calculation time/step < 0.1s for step<=10000
        surface_gradient_sum = np.zeros([2])
        for seq, loc in enumerate(target_points_abs_loc):
            gradient_history_part = calc_gradient_history_part_2d_optimized(loc, location_history, peclet_number, step, dt)
            gradient_local_part = calc_gradient_local_part_2d_optimized(loc, location_history, peclet_number, step, dt)
            chemical_gradient_sum = gradient_history_part + gradient_local_part
            surface_gradient = calc_tangential_gradient_part_numba(structure_ref_config, seq, chemical_gradient_sum)
            surface_gradient_sum += surface_gradient

        chem_gradient = 2 * np.pi/np.shape(structure_ref_config)[0] * surface_gradient_sum

    elif step >= 2 and acceleration == "parallel":
        # Accelerated with Python built-in parallel computing method apply_async
        # Parallel Pools: target points (discretization number of the structure)
        # Potentially useful for EXTRA-long simulation and fine grid of the structure (or even 3D cases)
        # Averaged calculation time/step ~= 0.1s
        surface_gradient_sum = np.zeros([2])
        target_points_num = np.arange(np.shape(structure_ref_config)[0])
        #num_cores = int(mp.cpu_count())
        process_pool = mp.Pool(num_core)
        #start_time = time.time()
        points_gradient_results = [process_pool.apply_async(calc_gradient_2d_parallel,
                                   args=(seq, target_points_abs_loc, location_history, structure_ref_config, peclet_number, step, dt))
                                   for seq in target_points_num]
        process_pool.close()
        process_pool.join()
        points_gradient = [segment.get() for segment in points_gradient_results]
        #elapsed_time = time.time() - start_time
        #print('parallel process time = ', elapsed_time)
        surface_gradient_sum += np.sum(points_gradient, 0)
        chem_gradient = 2 * np.pi/np.shape(structure_ref_config)[0] * surface_gradient_sum

    return chem_gradient


# --- Main ChemPhysics Calculation Function (for 2D multiple bodies) ---
def calc_surface_gradient_circle_numba_optimized(body_i, all_bodies, peclet_number, structure_ref_config, dt, *args, **kwargs):
    """
    Calculates the total surface gradient on a target body (body_i)
    by summing the contributions from ALL source bodies (all_bodies).

    Args:
        body_i: The target body for which to calculate the force.
        all_bodies: A list of all body objects in the simulation (the sources).
        peclet_number (Pe): The Peclet number.
        structure_ref_config: The vertex points defining the body's shape.
        dt: The time step size.
    """
    step = kwargs.get('step')

    # At step 0, there is no displacement and thus no gradient.
    if step == 0 and step >= 0:
        return np.zeros(2)

    target_points_abs_loc = structure_ref_config + body_i.location
    num_surface_points = np.shape(structure_ref_config)[0]
    total_surface_gradient_sum = np.zeros(2)

    # Loop over each surface point 's' on the target body_i
    for s_idx, s_loc in enumerate(target_points_abs_loc):

        chemical_gradient_at_point_s = np.zeros(2)

        # === CORE MULTI-BODY LOGIC ===
        # For this point 's', iterate through ALL source bodies (body_j)
        # to sum up the chemical field they generate at this location.
        for body_j in all_bodies:
            # At step 1, there is no "history" part. The gradient is only
            # determined by the "local" effect from the particle's state at step 0.
            if step == 1:
                # Call ONLY the local part calculation function.
                # It calculates the effect from the interval t=0 to t=1.
                grad_local_from_j = calc_gradient_local_part_2d_optimized(s_loc, body_j.location_history, peclet_number,
                                                                          step, dt)

                chemical_gradient_at_point_s += grad_local_from_j

            # For all subsequent steps (step > 1), we include both history and local parts.
            else:  # step > 1
                # The history part integrates from t=0 to t=(step-1).
                grad_hist_from_j = calc_gradient_history_part_2d_optimized(s_loc, body_j.location_history,
                                                                           peclet_number, step, dt)
                # The local part handles the most recent interval from t=(step-1) to t=step.
                grad_local_from_j = calc_gradient_local_part_2d_optimized(s_loc, body_j.location_history, peclet_number,
                                                                          step, dt)

                chemical_gradient_at_point_s += grad_hist_from_j + grad_local_from_j

        tangential_gradient_at_s = calc_tangential_gradient_part_numba(structure_ref_config, s_idx,
                                                                       chemical_gradient_at_point_s)

        total_surface_gradient_sum += tangential_gradient_at_s

    # Approximate the surface integral
    final_chem_force = (2 * np.pi / num_surface_points) * total_surface_gradient_sum
    return final_chem_force


def calc_surface_gradient_sphere_old(body, peclet_number, structure_ref_config, dt, *args, **kwargs):
    chem_gradient = np.zeros([3])
    step = kwargs.get('step')
    location = body.location
    location_history = body.location_history
    target_points_abs_loc = structure_ref_config + location
    acceleration = kwargs.get('acceleration')
    num_core = kwargs.get('core')

    # First guess
    if step == 0 and step >= 0:
        chem_gradient = np.zeros([3])
    # After initialized: second step only have
    elif step == 1 and step >= 0:
        surface_gradient_history_sum = np.zeros([3])
        for seq, loc in enumerate(target_points_abs_loc):
            gradient_history_part = calc_gradient_history_part_3d(loc, location_history, peclet_number, step, dt)
            surface_gradient_history = calc_tangential_gradient_part(structure_ref_config, seq, gradient_history_part)
            surface_gradient_history_sum += surface_gradient_history

        chem_gradient = 4 * np.pi/np.shape(structure_ref_config)[0] * surface_gradient_history_sum

    elif step >= 2 and acceleration == "numba":
        # Numba accelerated
        # Averaged calculation time/step < 0.1s for step<=10000
        surface_gradient_sum = np.zeros([3])
        for seq, loc in enumerate(target_points_abs_loc):
            gradient_history_part = calc_gradient_history_part_3d(loc, location_history, peclet_number, step, dt)
            # gradient_local_part = calc_gradient_local_part_3d(loc, location_history, peclet_number, step, dt)
            chemical_gradient_sum = gradient_history_part  # + gradient_local_part
            surface_gradient = calc_tangential_gradient_part(structure_ref_config, seq, chemical_gradient_sum)
            surface_gradient_sum += surface_gradient

        chem_gradient = 2 * np.pi/np.shape(structure_ref_config)[0] * surface_gradient_sum

    elif step >= 2 and acceleration == "parallel":
        # Accelerated with Python built-in parallel computing method apply_async
        # Parallel Pools: target points (discretization number of the structure)
        # Potentially useful for EXTRA-long simulation and fine grid of the structure
        # Averaged calculation time/step ~= 0.1s ?
        surface_gradient_sum = np.zeros([3])
        target_points_num = np.arange(np.shape(structure_ref_config)[0])
        #num_cores = int(mp.cpu_count())
        process_pool = mp.Pool(num_core)
        #start_time = time.time()
        points_gradient_results = [process_pool.apply_async(calc_gradient_2d_parallel,
                                   args=(seq, target_points_abs_loc, location_history, structure_ref_config, peclet_number, step, dt))
                                   for seq in target_points_num]
        process_pool.close()
        process_pool.join()
        points_gradient = [segment.get() for segment in points_gradient_results]
        #elapsed_time = time.time() - start_time
        #print('parallel process time = ', elapsed_time)
        surface_gradient_sum += np.sum(points_gradient, 0)
        chem_gradient = 4 * np.pi/np.shape(structure_ref_config)[0] * surface_gradient_sum

    return chem_gradient


def rotation_matrix_2d(theta):
    """
    Return the rotation matrix representing rotation
    by given an angle of rotation theta,
    which represents a rotation clockwise about the vector phi of magnitude phi.
    """
    return np.array([[np.cos(theta), -np.sin(theta)],
                    [np.sin(theta), np.cos(theta)]])


def rotation_matrix_3d(theta):
    """
    A 3D rotation encoded by an angle-axis representation as angle * axis
    see Rodrigues formula
    """
    theta_norm = np.linalg.norm(theta)
    theta_transform = np.array([[0, -theta[2], theta[1]],
                                [theta[2], 0, -theta[0]],
                                [-theta[1], theta[0], 0]])
    rotation_matrix = np.eye(3) + np.sin(theta_norm)/theta_norm * theta_transform \
                      + (1 - np.cos(theta_norm))/np.power(theta_norm, 2) * np.dot(theta_transform, theta_transform)

    return rotation_matrix


def sph2cart_field_3d(div, sph_vector):
    # sph2cart_field convert spherical coordinate to cartesian coordinate
    r, theta, phi = sph_vector[0], sph_vector[1], sph_vector[2]
    transform = np.array([[np.sin(theta) * np.cos(phi), np.cos(theta) * np.cos(phi), -np.sin(phi)],
                          [np.sin(theta) * np.cos(phi), np.cos(theta) * np.sin(phi), np.cos(phi)],
                          [np.cos(theta), -np.sin(theta), 0]])
    cart_vector = r * np.dot(transform, np.transpose(div))
    return np.transpose(cart_vector)


def cart2sph_vector_3d(cart_vector):
    # cart2sph_vector_3d convert cartesian coordinate to spherical coordinate
    r = np.linalg.norm(cart_vector)
    theta = np.arccos(cart_vector[2] / r)
    phi = np.arctan2(cart_vector[1], cart_vector[0])

    sph_vector = np.array([r, theta, phi])
    return sph_vector


# ==============================================================================
# 3D Function (Using Analytical Formula for Local Part)
# ==============================================================================

def get_chem_grad_3D(body, peclet_number, dt, *args, **kwargs):
    """
    Calculates the 3D chemical gradient.
    - History part: Sum over discrete time steps.
    - Local part: Uses the analytical formula from the paper (Eq. 17).
    """
    surface_nodes = body.get_surface_nodes()
    total_grad_C_on_nodes = np.zeros_like(surface_nodes)
    D = peclet_number ** -1
    step = kwargs.get('step')
    t_now = step * dt

    # History Part (C_H): Sum over discrete steps from t=0 to t_now - dt
    if step > 1:
        for i in range(step - 1):
            t_prime = i * dt
            tau = t_now - t_prime
            pos_prime = body.location_history[i]

            for j, r_s in enumerate(surface_nodes):
                r_vec = r_s - pos_prime
                r_sq = np.dot(r_vec, r_vec)
                if tau <= 1e-12: continue
                G = (4 * np.pi * D * tau) ** (-1.5) * np.exp(-r_sq / (4 * D * tau))
                grad_G = -r_vec / (2 * D * tau) * G
                total_grad_C_on_nodes[j] += grad_G * dt

    # Local Part (C_L): Analytical solution for the integral over the last time step
    if step > 0:
        pos_t_now = body.location
        pos_t_before = body.location_history[-1]
        v_k = (pos_t_now - pos_t_before) / dt
        V_sq = np.dot(v_k, v_k)

        for j, r_s in enumerate(surface_nodes):
            r_sk_vec = r_s - pos_t_now
            R = np.linalg.norm(r_sk_vec)

            if R < 1e-12: continue

            # Implementation of analytical solution from Eq. 15, 16, 17
            r_sk_hat = r_sk_vec / R
            v_dot_r_hat = np.dot(v_k, r_sk_hat)
            v_parallel_vec = v_dot_r_hat * r_sk_hat
            v_perp_vec = v_k - v_parallel_vec
            v_perp_sq = np.dot(v_perp_vec, v_perp_vec)

            arg_erfc = -v_dot_r_hat / (2 * np.sqrt(D / dt))

            exp_term_F = np.exp(-V_sq * dt / (4 * D))
            exp_term_R = np.exp(-R ** 2 / (4 * D * dt))
            F = (1 / R ** 3) * (exp_term_F * erfc(arg_erfc) - exp_term_R)

            exp_term_G = np.exp(-v_perp_sq * dt / (4 * D))
            G = (1 / R) * np.sqrt(dt / (np.pi * D)) * exp_term_G

            grad_C_local = (1 / (4 * np.pi * D)) * (F * r_sk_vec + G * v_k)
            total_grad_C_on_nodes[j] += grad_C_local

    return total_grad_C_on_nodes


# ==============================================================================
# 3D Functions (Numba-accelerated version with Analytical Local Part)
# ==============================================================================

@njit(parallel=True, fastmath=True)
def _history_part_3d_numba(n_nodes, surface_nodes, n_history_steps, pos_history_arr, t_now, dt, D):
    """
    Numba-jitted function for the history part of the gradient.
    NOTE: This function's formula was correct. It calculates the gradient of the
    concentration's Green's function (grad_G), not the concentration itself.
    """
    total_grad_C = np.zeros((n_nodes, 3))
    for j in prange(n_nodes):
        r_s = surface_nodes[j]
        grad_sum = np.zeros(3)
        for i in range(n_history_steps):
            t_prime = i * dt
            tau = t_now - t_prime
            pos_prime = pos_history_arr[i]
            r_vec = r_s - pos_prime
            r_sq = r_vec[0] ** 2 + r_vec[1] ** 2 + r_vec[2] ** 2
            if tau <= 1e-12: continue
            # G is the Green's function for concentration
            const = (4 * np.pi * D * tau)
            G = np.power(const, -1.5) * np.exp(-r_sq / const)
            # grad_G is the gradient of G, which is what we need
            grad_G = -r_vec / (2 * D * tau) * G
            grad_sum += grad_G
        # Integration using Euler method
        total_grad_C[j] = grad_sum * dt
    return total_grad_C


@njit(parallel=True, fastmath=True)
def _local_part_3d_numba_numerical(n_nodes, surface_nodes, pos_t_now, pos_t_before, Pe, dt):
    """
    Numba-jitted function for the local part of the gradient, using numerical integration (Midpoint Rule).
    This replaces the incorrect analytical formula.
    """
    total_grad_C_local = np.zeros((n_nodes, 3))
    v_k = (pos_t_now - pos_t_before) / dt

    # Constant factor Eq. E2 for d=3
    # -2 / (pi^(3/2)) * (Pe/4)^(5/2) = -Pe^(2.5) / (16 * pi^(1.5))
    const_factor = -np.power(Pe, 2.5) / (16.0 * np.power(np.pi, 1.5))

    # Midpoint rule for the local integral over [t_now-dt, t_now]
    # Time point is t_now - dt/2, so tau = t_now - t' = dt/2
    tau_mid = 0.5 * dt
    # Position at midpoint time, assuming constant velocity in the last step
    r_prime_mid = pos_t_before + v_k * (0.5 * dt)

    for j in prange(n_nodes):
        r_s = surface_nodes[j]
        r_vec = r_s - r_prime_mid
        r_sq = r_vec[0] ** 2 + r_vec[1] ** 2 + r_vec[2] ** 2

        # Integrand from Eq. E2 evaluated at the midpoint
        integrand_val = const_factor * (r_vec / np.power(tau_mid, 2.5)) * np.exp(-r_sq * Pe / (4.0 * tau_mid))

        # The integral is approximately the integrand at midpoint * interval width (dt)
        total_grad_C_local[j] = integrand_val * dt

    return total_grad_C_local


def get_chem_grad_3D_numba(body, peclet_number, dt, *args, **kwargs):
    """Numba-accelerated version using numerical integration for the local part."""
    surface_nodes = body.get_surface_nodes(body.location, body.omega_axis_orientation)
    total_grad_C_on_nodes = np.zeros_like(surface_nodes)
    D = 1.0 / peclet_number
    step = kwargs.get('step')
    t_now = step * dt
    n_nodes = body.n_nodes

    # History Part
    if step > 1:
        # History is from t=0 to t=(step-2)*dt
        n_history_steps = step - 1
        pos_history_arr = np.array(body.location_history[:n_history_steps])
        grad_C_history = _history_part_3d_numba(n_nodes, surface_nodes, n_history_steps, pos_history_arr, t_now, dt, D)
        total_grad_C_on_nodes += grad_C_history

    # Local Part (Numerical)
    if step > 0:
        pos_t_now = body.location
        # The last position in history is at step-1
        pos_t_before = body.location_history[step - 1]
        grad_C_local = _local_part_3d_numba_numerical(n_nodes, surface_nodes, pos_t_now, pos_t_before, peclet_number,
                                                      dt)
        total_grad_C_on_nodes += grad_C_local

    return total_grad_C_on_nodes


def calc_tangential_grad_3D(body, peclet_number, dt, *args, **kwargs):
    """
    This is the main public function to get the tangential chemical gradient on the body surface.
    It acts as a wrapper, calling the appropriate backend (SciPy or Numba) to get the
    full 3D gradient, and then projects it to the surface to get the tangential component.

    Args:
        body (Body3D): The 3D body object.
        sim_config (object): The simulation configuration object.
        use_numba (bool): Flag to select the calculation backend. True for Numba, False for SciPy.

    Returns:
        np.ndarray: An array of tangential gradient vectors for each surface node.
    """
    # Step 1: Get the full 3D gradient using the selected backend
    acceleration = kwargs.get('acceleration')
    if acceleration == "numba":
        grad_C_on_nodes = get_chem_grad_3D_numba(body, peclet_number, dt, *args, **kwargs)
    else:
        grad_C_on_nodes = get_chem_grad_3D(body, peclet_number, dt, *args, **kwargs)

    # Step 2: Project the full gradient onto the tangential plane (∇sC)
    # =================================================================
    # Get the current world coordinates of the surface nodes
    nodes_world = body.get_surface_nodes(body.location, body.omega_axis_orientation)

    # Calculate normal vectors (from center to each node)
    r_vectors = nodes_world - body.location
    r_vectors_norm = np.linalg.norm(r_vectors, axis=1, keepdims=True)
    # Avoid division by zero for any node at the center (theoretically shouldn't happen)
    r_vectors_norm[r_vectors_norm < 1e-12] = 1.0
    normals = r_vectors / r_vectors_norm

    # Project ∇C onto the normal vectors to get the normal component
    grad_dot_norm = np.sum(grad_C_on_nodes * normals, axis=1, keepdims=True)
    grad_normal_component = grad_dot_norm * normals

    # The tangential gradient is the full gradient minus its normal component
    grad_tangential_on_nodes = grad_C_on_nodes - grad_normal_component

    # Step 3: Integrate over the surface and calculate the final force (Fc)
    # ======================================================================
    # The surface integral is approximated by taking the mean of the tangential
    # gradient vectors over all surface nodes.
    mean_tangential_grad = np.mean(grad_tangential_on_nodes, axis=0)

    return mean_tangential_grad


# ==============================================================================
# HELPER FUNCTIONS for cases: dry point case and spherical particle case
# ==============================================================================

@njit(fastmath=True)
def quaternion_to_rotation_matrix(q):
    """
    Quaternion (w, x, y, z) transform into a 3x3 rotation matrix.
    Warning：for numba capacity，the input sequence is [x, y, z, w].
    """
    x, y, z, w = q[0], q[1], q[2], q[3]

    # Pre-calculate squared terms
    x2, y2, z2 = x * x, y * y, z * z

    # Pre-calculate products
    xy, xz, yz = x * y, x * z, y * z
    wx, wy, wz = w * x, w * y, w * z

    # 构建旋转矩阵
    R = np.empty((3, 3))
    R[0, 0] = 1.0 - 2.0 * (y2 + z2)
    R[0, 1] = 2.0 * (xy - wz)
    R[0, 2] = 2.0 * (xz + wy)

    R[1, 0] = 2.0 * (xy + wz)
    R[1, 1] = 1.0 - 2.0 * (x2 + z2)
    R[1, 2] = 2.0 * (yz - wx)

    R[2, 0] = 2.0 * (xz - wy)
    R[2, 1] = 2.0 * (yz + wx)
    R[2, 2] = 1.0 - 2.0 * (x2 + y2)

    return R


# ==============================================================================
# Case 1：POINT CASE (Structures only translate NO rotation)
# ==============================================================================

@njit(parallel=True, fastmath=True)
def _history_part_point_3d_numba(n_nodes, surface_nodes, n_history_steps, pos_history_arr, t_now, dt, D):
    """
    Calculate history part for point particle case.
    sigma=1 distribute uniformly.
    """
    total_grad_C = np.zeros((n_nodes, 3))
    for j in prange(n_nodes):
        # r_s is position in laboratory frame
        r_s = surface_nodes[j]
        grad_sum = np.zeros(3)
        for i in range(n_history_steps):
            t_prime = i * dt
            tau = t_now - t_prime
            pos_prime_center = pos_history_arr[i]
            r_vec = r_s - pos_prime_center
            r_sq = r_vec[0] ** 2 + r_vec[1] ** 2 + r_vec[2] ** 2
            if tau <= 1e-12: continue
            const = (4 * np.pi * D * tau)
            G = np.power(const, -1.5) * np.exp(-r_sq / const)
            grad_G = -r_vec / (2 * D * tau) * G
            grad_sum += grad_G
        # Integration using Euler method
        total_grad_C[j] = grad_sum * dt
    return total_grad_C


@njit(parallel=True, fastmath=True)
def _local_part_point_3d_numba(n_nodes, surface_nodes, pos_t_now, pos_t_before, Pe, dt):
    """
    Calculate local part for point particle case.
    """
    total_grad_C_local = np.zeros((n_nodes, 3))
    v_k = (pos_t_now - pos_t_before) / dt

    const_factor = -np.power(Pe, 2.5) / (16.0 * np.power(np.pi, 1.5))
    tau_mid = 0.5 * dt
    r_prime_mid = pos_t_before + v_k * (0.5 * dt)

    for j in prange(n_nodes):
        r_s = surface_nodes[j]
        r_vec = r_s - r_prime_mid
        r_sq = r_vec[0] ** 2 + r_vec[1] ** 2 + r_vec[2] ** 2

        integrand_val = const_factor * (r_vec / np.power(tau_mid, 2.5)) * np.exp(-r_sq * Pe / (4.0 * tau_mid))
        total_grad_C_local[j] = integrand_val * dt

    return total_grad_C_local


def history_local_compose_3d_point(body, peclet_number, dt, *args, **kwargs):
    """
    Point Case: combining both local and history parts
    Only return the chemical force.
    """
    surface_nodes = body.get_surface_nodes(body.location, body.omega_axis_orientation, body.is_janus)
    total_grad_C_on_nodes = np.zeros_like(surface_nodes)
    D = 1.0 / peclet_number
    step = kwargs.get('step')
    t_now = step * dt
    n_nodes = body.n_nodes

    # History Part
    if step > 1:
        n_history_steps = step - 1
        pos_history_arr = np.array(body.location_history[:n_history_steps])
        grad_C_history = _history_part_point_3d_numba(
            n_nodes, surface_nodes, n_history_steps, pos_history_arr, t_now, dt, D)
        total_grad_C_on_nodes += grad_C_history

    # Local Part
    if step > 0:
        pos_t_now = body.location
        # The last position in history is at step-1
        pos_t_before = body.location_history[step - 1]
        grad_C_local = _local_part_point_3d_numba(
            n_nodes, surface_nodes, pos_t_now, pos_t_before, peclet_number, dt)
        total_grad_C_on_nodes += grad_C_local

    # 计算切向梯度和最终的力
    normals = surface_nodes - body.location
    normals /= np.linalg.norm(normals, axis=1)[:, np.newaxis]
    tangential_grad = total_grad_C_on_nodes - (np.sum(total_grad_C_on_nodes * normals, axis=1))[:, np.newaxis] * normals
    mean_tangential_grad_force = np.mean(tangential_grad, axis=0)

    return mean_tangential_grad_force


# ==============================================================================================================
# Case 2：Spherical particle with surface chemical substance distribution (Janus particle, with body frame rotation)
# ==============================================================================================================
@njit(parallel=True, fastmath=True)
def _precompute_rotated_nodes(n_history_steps, n_nodes, pos_history, orient_history, nodes_body_frame):
    """
    Precompute rotated node positions for history time steps
    """
    rotated_nodes_history = np.empty((n_history_steps, n_nodes, 3))
    for i in prange(n_history_steps):
        pos_center_hist = pos_history[i]
        orient_hist = orient_history[i]
        rot_matrix_hist = quaternion_to_rotation_matrix(orient_hist)
        rotated_nodes_history[i, :, :] = pos_center_hist + np.dot(nodes_body_frame, rot_matrix_hist.T)
    return rotated_nodes_history


@njit(parallel=True, fastmath=True)
def _history_part_dist_3d_numba_optimized(n_nodes, rotated_nodes_history, sigma_dist, t_now, dt, D):
    """
    Distribution Case:
    Calculate history part contribution,
    with consideration of particle self rotation and distribution of chemical substance.
    Optimized with precompute_rotated_nodes.
    """
    total_grad_C = np.zeros((n_nodes, 3))
    n_history_steps = rotated_nodes_history.shape[0]

    sigma_dist_scaled = sigma_dist.reshape(-1, 1)

    for j in prange(n_nodes):
        # location of target node in body frame
        target_node_lab_hist_all_steps = rotated_nodes_history[:, j, :]
        grad_sum = np.zeros(3)

        for i in range(n_history_steps):
            target_node_lab_hist = target_node_lab_hist_all_steps[i]
            # Retrieve all source nodes positions in history
            source_nodes_lab_hist = rotated_nodes_history[i, :, :]

            # Vectorized computation: vector distance and distance square for all source nodes to target nodes
            r_vec_all_sources = target_node_lab_hist - source_nodes_lab_hist
            r_sq_all_sources = np.sum(r_vec_all_sources**2, axis=1)

            t_prime = i * dt
            tau = t_now - t_prime

            if tau > 1e-12:
                const = (4 * np.pi * D * tau)
                G_all_sources = np.power(const, -1.5) * np.exp(-r_sq_all_sources / const)
                grad_G_all_sources = (-r_vec_all_sources / (2 * D * tau)) * G_all_sources.reshape(-1, 1)

                # Vectorized computation: gradient multiplied by sigma value
                weighted_grads = grad_G_all_sources * sigma_dist_scaled
                grad_sum += np.sum(weighted_grads, axis=0)

        total_grad_C[j] = grad_sum * dt
    return total_grad_C


@njit(parallel=True, fastmath=True)
def _local_part_dist_3d_numba(n_nodes, pos_now, orient_now, pos_before, orient_before, nodes_body_frame, sigma_dist, Pe,
                              dt):
    """
    Distribution Case:
    Calculate history part contribution,
    using mid-point rule for handling body rotation.
    """
    total_grad_C_local = np.zeros((n_nodes, 3))
    D = 1.0 / Pe

    # mid-point for position and orientation
    pos_mid = (pos_now + pos_before) / 2.0
    # nlerp is a good approximation for small dt
    orient_mid = (orient_now + orient_before) / 2.0
    orient_mid /= np.linalg.norm(orient_mid)
    rot_matrix_mid = quaternion_to_rotation_matrix(orient_mid)

    for j in prange(n_nodes):
        target_node_body = nodes_body_frame[j]
        target_node_lab_mid = pos_mid + np.dot(rot_matrix_mid, target_node_body)

        grad_sum = np.zeros(3)
        for k in range(n_nodes):
            sigma_k = sigma_dist[k]
            if sigma_k == 0: continue

            source_node_body = nodes_body_frame[k]
            source_node_lab_mid = pos_mid + np.dot(rot_matrix_mid, source_node_body)

            tau = 0.5 * dt
            r_vec = target_node_lab_mid - source_node_lab_mid
            r_sq = r_vec[0] ** 2 + r_vec[1] ** 2 + r_vec[2] ** 2

            if r_sq < 1e-12: continue

            const = (4 * np.pi * D * tau)
            G = np.power(const, -1.5) * np.exp(-r_sq / const)
            grad_G = -r_vec / (2 * D * tau) * G
            grad_sum += grad_G * sigma_k

        total_grad_C_local[j] = grad_sum * dt

    return total_grad_C_local


@njit(parallel=True, fastmath=True)
def _local_part_dist_3d_numba_optimized(n_nodes, pos_now, orient_now, pos_before, orient_before, nodes_body_frame,
                                        sigma_dist, Pe, dt):
    """
    Local Part: vectorized version for optimized computation time
    """
    total_grad_C_local = np.zeros((n_nodes, 3))
    D = 1.0 / Pe

    # mid position and orientation for mid-point rule (CoM)
    pos_mid = (pos_now + pos_before) / 2.0
    orient_mid = (orient_now + orient_before) / 2.0
    orient_mid /= np.linalg.norm(orient_mid)
    rot_matrix_mid = quaternion_to_rotation_matrix(orient_mid)

    # mid position and orientation for all nodes
    nodes_lab_mid = pos_mid + np.dot(nodes_body_frame, rot_matrix_mid.T)

    # sigma_dist is a (n_nodes,) 1d array
    sigma_dist_scaled = sigma_dist.reshape(-1, 1)

    tau = 0.5 * dt
    const = (4 * np.pi * D * tau)

    for j in prange(n_nodes):
        target_node_lab_mid = nodes_lab_mid[j]

        # Vectorized computations: vec distance and distance square for all source nodes to target nodes.
        r_vecs = target_node_lab_mid - nodes_lab_mid
        r_sqs = np.sum(r_vecs ** 2, axis=1)

        # Calculate Greens' function and its gradient for all source nodes.
        # avoiding dividing by zero when it is too close.
        valid_indices = r_sqs > 1e-12

        G = np.zeros(n_nodes)
        G[valid_indices] = np.power(const, -1.5) * np.exp(-r_sqs[valid_indices] / const)

        grad_G = -r_vecs / (2 * D * tau) * G.reshape(-1, 1)

        # gradient multiplied by sigma value.
        weighted_grad_sum = np.sum(grad_G * sigma_dist_scaled, axis=0)
        total_grad_C_local[j] = weighted_grad_sum * dt

    return total_grad_C_local


def history_local_compose_3d_distribution(body, peclet_number, dt, *args, **kwargs):
    """
    Distribution case:
    Combining history and local parts and calculate chemical gradients for chemical force and torque.
    RETURN:
        chemical force, chemical torque
    """
    step = kwargs.get('step')
    t_now = step * dt
    n_nodes = body.n_nodes
    D = 1.0 / peclet_number

    nodes_body_frame = body.nodes_body_frame
    sigma_dist = body.sigma_distribution

    total_grad_C_on_nodes = np.zeros((n_nodes, 3))

    # History Part
    if step > 1:
        n_history_steps = step - 1
        pos_history = np.array(body.location_history[:n_history_steps])
        orient_history = np.array(body.orientation_history[:n_history_steps])

        # Call helper function to precompute all nodes
        rotated_nodes_history = _precompute_rotated_nodes(n_history_steps, n_nodes, pos_history, orient_history,
                                                          nodes_body_frame)

        # Call optimized history part function
        grad_C_history = _history_part_dist_3d_numba_optimized(
            n_nodes, rotated_nodes_history, sigma_dist, t_now, dt, D)
        total_grad_C_on_nodes += grad_C_history

    # Local Part
    if step > 0:
        pos_now = body.location
        orient_now = body.orientation
        pos_before = body.location_history[step - 1]
        orient_before = body.orientation_history[step - 1]

        # Call optimized local part function
        grad_C_local = _local_part_dist_3d_numba_optimized(
            n_nodes, pos_now, orient_now, pos_before, orient_before, nodes_body_frame, sigma_dist, peclet_number, dt)
        total_grad_C_on_nodes += grad_C_local

    # Calculate chemical force and torque
    surface_nodes = body.get_surface_nodes()
    normals = surface_nodes - body.location
    normals /= np.linalg.norm(normals, axis=1)[:, np.newaxis]

    tangential_grad = total_grad_C_on_nodes - (np.sum(total_grad_C_on_nodes * normals, axis=1))[:, np.newaxis] * normals

    # mobility should be multiplied in the integrator_3D.py
    mean_tangential_grad_force = np.mean(tangential_grad, axis=0)

    # toque is integration of r_j x F_j
    # F_j = -mobility * tangential_grad[j]
    # Torque_j = r_j x F_j = (nodes_body_frame[j]) x (-mobility * tangential_grad[j])
    # mobility should be multiplied in the integrator_3D.py
    r_vectors = nodes_body_frame
    torque_vectors = np.cross(r_vectors, tangential_grad)
    mean_torque = np.mean(torque_vectors, axis=0)

    return mean_tangential_grad_force, mean_torque
