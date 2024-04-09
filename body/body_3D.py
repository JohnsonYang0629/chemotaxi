import numpy as np
import copy
import sys

class Body_3D(object):
  '''
  Small class to handle a single body in 3D domain.
  '''
  def __init__(self, location, orientation):
    '''
    Constructor. Take arguments like ...
    '''
    # Location as np.array.shape = 2
    self.location = np.copy(location)
    self.location_new = np.copy(location)
    self.location_old = np.copy(location)
    # Orientation as Quaternion
    self.orientation = copy.copy(orientation)
    self.orientation_new = copy.copy(orientation)
    self.orientation_old = copy.copy(orientation)
    # Reference configuration. Coordinates of droplet for quaternion [1, 0, 0, 0]
    # and location = np.array[0, 0, 0]) as a np.array.shape = (1, 3)
    # self.rotation_matrix = None
    # Some default functions
    self.function_force = self.default_none
    self.function_torque = self.default_none
    self.prescribed_velocity = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
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