import numpy as np
import copy
import sys


class Body2D(object):
  '''
  Small class to handle a single body in 2D domain.
  '''
  def __init__(self, location, orientation, n_steps):
    '''
    Constructor. Take arguments like ...
    '''
    # Location as np.array.shape = 2
    self.location = np.copy(location)
    self.location_new = np.copy(location)
    self.location_old = np.copy(location)
    self.location_history = np.zeros([n_steps + 1, 2])
    # Orientation as Quaternion
    self.orientation = np.copy(orientation)
    self.orientation_new = np.copy(orientation)
    self.orientation_old = np.copy(orientation)
    # Reference configuration. Coordinates of droplet for polar coordinate theta [0]
    # and location = np.array[0, 0]) as a np.array.shape = (1, 2)
    # self.rotation_matrix = None
    # Some default functions
    self.function_force = self.default_none
    self.function_torque = self.default_none
    self.prescribed_velocity = np.array([0.0, 0.0, 0.0])
    self.chem_surface_gradient = np.array([0.0, 0.0])
    self.grad_chem_local_part = np.array([0.0, 0.0])
    self.grad_chem_history_part = np.array([0.0, 0.0])
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
