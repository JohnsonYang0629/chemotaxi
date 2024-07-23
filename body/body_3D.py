import numpy as np
import copy
import sys


class Body3D(object):
  '''
  Small class to handle a single body in 3D domain.
  '''
  def __init__(self, location, v_orientation, omega_orientation, n_steps):
    '''
    Constructor. Take arguments like ...
    '''
    # Location as np.array.shape = 3
    self.location = np.copy(location)
    self.location_new = np.copy(location)
    self.location_old = np.copy(location)
    self.location_history = np.zeros([n_steps + 1, 3])
    # Orientation using Rodrigues formula
    self.v_orientation = copy.copy(v_orientation)
    self.v_orientation_new = copy.copy(v_orientation)
    self.v_orientation_old = copy.copy(v_orientation)
    self.omega_orientation = copy.copy(omega_orientation)
    self.omega_orientation_new = copy.copy(omega_orientation)
    self.omega_orientation_old = copy.copy(omega_orientation)
    # Reference configuration. Coordinates of droplet for quaternion [1, 0, 0, 0]
    # and location = np.array[0, 0, 0]) as a np.array.shape = (1, 3)
    # Some default functions
    self.function_force = self.default_none
    self.function_torque = self.default_none
    self.prescribed_velocity = np.array([0.0, 0.0, 0.0, 0.0])
    self.chem_surface_gradient = np.array([0.0, 0.0, 0.0])
    self.ID = None

  def calc_prescribed_velocity(self):
      '''
      Return the body prescribed velocity.
      '''
      return self.prescribed_velocity

  def calc_force(self):
    '''
    Return the force on the body.
    '''
    return self.function_force()

  def calc_torque(self):
    '''
    Return the torque on the body.
    '''
    return self.function_torque()

  def default_none(self, *args, **kwargs):
    return None