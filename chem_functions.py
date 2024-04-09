import numpy as np
import sys


def calc_surface_gradient_circle(body, peclet_number, structure_ref_config, *args, **kwargs):
    chem_gradient = np.zeros([1, 2])
    step = kwargs.get('step')
    if step == 0:
        chem_gradient = np.zeros([1, 2])
    else:
        chem_gradient = np.zeros([1, 2])

    return chem_gradient


def rotation_matrix_2d(theta):
    """
    Return the rotation matrix representing rotation
    by given an angle of rotation theta,
    which represents a rotation clockwise about the vector phi of magnitude phi.
    """
    return np.array([[np.cos(theta), -np.sin(theta)],
                    [np.sin(theta), np.cos(theta)]])
