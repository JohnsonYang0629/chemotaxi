import numpy as np
import scipy.special as sc
import sys


def calc_gradient_history_part_2d(target_position, history_path_location, peclet_number, step, dt):
    gradient_segment = np.zeros([step, 2])
    for s in range(step):
        position_difference = target_position - history_path_location[s, :]
        position_difference_norm_square = np.power(np.linalg.norm(position_difference), 2)
        gradient_segment[s, :] = 2 * np.power(peclet_number/4, 2)/(np.pi * np.power((step - s) * dt, 2)) \
                                   * np.exp(-peclet_number/(4 * (step - s) * dt) * position_difference_norm_square)\
                                   * position_difference
    # Use trapezoidal rule to sum up segments
    gradient_history_part = 1/2 * dt * np.sum(gradient_segment)
    return gradient_history_part


def calc_gradient_local_part_2d(target_position, history_path_location, peclet_number, step, dt):
    position_difference = target_position - history_path_location[step - 1, :]
    position_difference_norm_square = np.power(np.linalg.norm(position_difference), 2)
    rho = peclet_number * position_difference_norm_square/(4 * dt)
    gradient_local_part = 2 * peclet_number/(4 * np.pi * np.power(position_difference_norm_square, 2)) \
                            * np.exp(-rho) * (np.power(rho, 2) + 2 * rho + 2) \
                            * position_difference
    return gradient_local_part


def calc_tangential_gradient_part(structure_ref_config, seq, gradient):
    pos = structure_ref_config[seq, :]
    pos_norm = pos/np.linalg.norm(pos)
    # Given that discretized points are on a circle
    pos_tangent_vec = np.array([pos_norm[1], -pos_norm[0]])
    tangential_gradient_part = np.dot(gradient, pos_tangent_vec) * pos_tangent_vec
    return tangential_gradient_part


def calc_surface_gradient_circle(body, peclet_number, structure_ref_config, dt, *args, **kwargs):
    chem_gradient = np.zeros([2])
    step = kwargs.get('step')
    location = body.location
    location_history = body.location_history
    target_points_abs_loc = structure_ref_config + location

    # first guess
    if step == 0 and step >= 0:
        chem_gradient = np.zeros([2])
    # after initialized: second step only have
    elif step == 1 and step >= 0:
        surface_gradient_history_sum = np.zeros([2])
        for seq, loc in enumerate(target_points_abs_loc):
            gradient_history_part = calc_gradient_history_part_2d(loc, location_history, peclet_number, step, dt)
            surface_gradient_history = calc_tangential_gradient_part(structure_ref_config, seq, gradient_history_part)
            surface_gradient_history_sum += surface_gradient_history

        chem_gradient = surface_gradient_history_sum
    elif step >= 2 and step >= 0:
        surface_gradient_sum = np.zeros([2])
        for seq, loc in enumerate(target_points_abs_loc):
            gradient_history_part = calc_gradient_history_part_2d(loc, location_history, peclet_number, step, dt)
            gradient_local_part = calc_gradient_local_part_2d(loc, location_history, peclet_number, step, dt)
            surface_gradient = calc_tangential_gradient_part(structure_ref_config, seq, (gradient_history_part + gradient_local_part))
            surface_gradient_sum += surface_gradient

        chem_gradient = surface_gradient_sum
    return chem_gradient


def rotation_matrix_2d(theta):
    """
    Return the rotation matrix representing rotation
    by given an angle of rotation theta,
    which represents a rotation clockwise about the vector phi of magnitude phi.
    """
    return np.array([[np.cos(theta), -np.sin(theta)],
                    [np.sin(theta), np.cos(theta)]])
