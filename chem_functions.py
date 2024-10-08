import numpy as np
import scipy.special as sc
import sys
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
    surface_gradient_local = calc_tangential_gradient_part(structure_ref_config, seq, gradient_local_part)
    surface_gradient_history = calc_tangential_gradient_part(structure_ref_config, seq, gradient_history_part)
    return [surface_gradient, surface_gradient_local, surface_gradient_history]


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
    #gradient_local_part = 2 * peclet_number/(4 * np.pi * np.power(position_difference_norm_square, 2)) \
    #                        * np.exp(-rho) * (np.power(rho, 2) + 2 * rho + 2) \
    #                        * position_difference
    gradient_local_part = 2 * peclet_number / (4 * np.pi * position_difference_norm_square) \
                            * np.exp(-rho) * position_difference
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
    # 2D case for planar particle (a circle)
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
        chem_local_part = np.zeros([2])
        chem_history_part = np.zeros([2])
    # After initialized: second step only have
    elif step == 1 and step >= 0:
        surface_gradient_history_sum = np.zeros([2])
        for seq, loc in enumerate(target_points_abs_loc):
            gradient_history_part = calc_gradient_history_part_2d(loc, location_history, peclet_number, step, dt)
            surface_gradient_history = calc_tangential_gradient_part(structure_ref_config, seq, gradient_history_part)
            surface_gradient_history_sum += surface_gradient_history

        chem_gradient = 2 * np.pi / np.shape(structure_ref_config)[0] * surface_gradient_history_sum
        chem_local_part = np.zeros([2])
        chem_history_part = 2 * np.pi / np.shape(structure_ref_config)[0] * surface_gradient_history_sum

    elif step >= 2 and acceleration == "numba":
        # Numba accelerated
        # Averaged calculation time/step < 0.1s for step<=10000
        surface_gradient_sum = np.zeros([2])
        surface_gradient_local_sum = np.zeros([2])
        surface_gradient_history_sum = np.zeros([2])
        for seq, loc in enumerate(target_points_abs_loc):
            gradient_history_part = calc_gradient_history_part_2d(loc, location_history, peclet_number, step, dt)
            gradient_local_part = calc_gradient_local_part_2d(loc, location_history, peclet_number, step, dt)
            chemical_gradient_sum = gradient_history_part + gradient_local_part
            surface_gradient = calc_tangential_gradient_part(structure_ref_config, seq, chemical_gradient_sum)
            surface_gradient_sum += surface_gradient
            surface_gradient_local_sum += calc_tangential_gradient_part(structure_ref_config, seq, gradient_local_part)
            surface_gradient_history_sum += calc_tangential_gradient_part(structure_ref_config, seq, gradient_history_part)

        chem_gradient = 2 * np.pi / np.shape(structure_ref_config)[0] * surface_gradient_sum
        chem_local_part = 2 * np.pi / np.shape(structure_ref_config)[0] * surface_gradient_local_sum
        chem_history_part = 2 * np.pi / np.shape(structure_ref_config)[0] * surface_gradient_history_sum

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
        points_gradient = np.sum(np.array(points_gradient), 0)
        chem_local_part = points_gradient[1]
        chem_history_part = points_gradient[2]
        #elapsed_time = time.time() - start_time
        #print('parallel process time = ', elapsed_time)
        surface_gradient_sum += points_gradient[0]
        chem_gradient = 2 * np.pi / np.shape(structure_ref_config)[0] * surface_gradient_sum

    return chem_gradient, chem_local_part, chem_history_part


def calc_surface_gradient_sphere(body, peclet_number, structure_ref_config, dt, *args, **kwargs):
    # 3D case for spherical particle
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


