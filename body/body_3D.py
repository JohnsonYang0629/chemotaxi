import numpy as np
import copy
import sys
from read_input import read_vertex_file
from integrator.quaternion import Quaternion


class Body3D(object):
  '''
  Small class to handle a single body in 3D domain.
  '''
  def __init__(self, location, omega_axis_orientation, structure_ref_config, n_steps):
    '''
    Constructor. Take arguments like ...
    '''
    # Location as np.array.shape = 3
    self.location = np.copy(location)
    self.location_new = np.copy(location)
    self.location_old = np.copy(location)
    self.location_history = np.zeros([n_steps + 1, 3])
    # Orientation as Quaternion
    self.omega_axis_orientation = copy.copy(omega_axis_orientation)
    self.omega_axis_orientation_new = copy.copy(omega_axis_orientation)
    self.omega_axis_orientation_old = copy.copy(omega_axis_orientation)
    # v0 axis as in 3D vector
    self.v0_axis = np.array([1.0, 0.0, 0.0])
    self.v0_axis_new = np.array([1.0, 0.0, 0.0])
    self.v0_axis_old = np.array([1.0, 0.0, 0.0])

    # omega axis as in 3D vector
    self.omega_axis = np.array([0.0, 0.0, 1.0])
    self.omega_axis_new = np.array([0.0, 0.0, 1.0])
    self.omega_axis_old = np.array([0.0, 0.0, 1.0])
    # Reference configuration. Coordinates of droplet for quaternion [1, 0, 0, 0]
    # and location = np.array[0, 0, 0]) as a np.array.shape = (1, 3)
    self.orientation = np.zeros((1, 4))
    self.orientation_history = np.zeros((n_steps + 1, 4))

    # Load surface node positions
    # These nodes are defined in the body's own reference frame.
    self.nodes_body_frame = np.copy(structure_ref_config)
    self.n_nodes = len(self.nodes_body_frame)

    # Particle surface chemical substance distribution
    self.sigma_distribution = np.ones(self.n_nodes) / self.n_nodes
    self.is_janus = False

    self.function_force = self.default_none
    self.function_torque = self.default_none
    self.prescribed_velocity = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    self.chem_surface_gradient = np.array([0.0, 0.0, 0.0])
    self.chem_torque_gradient = np.array([0.0, 0.0, 0.0])
    self.ID = None

  def set_sigma_distribution(self, sigma_values):
    """
    Set and normalize sigma disribution
    Args:
        sigma_values (list or np.ndarray): each node's corresponding sigma value.
    """
    # 1. check the length of sigma value is identical to the nodes.
    if sigma_values is None or len(sigma_values) != self.n_nodes:
      print("Warning: Input sigma_values invalid，default uniform distribution will be used.")
      # if costume distribution invalid, we use non-Janus particle.
      self.sigma_distribution = np.ones(self.n_nodes) / self.n_nodes
      self.is_janus = False
      return

    # 2. normalize with the sum of total sigma
    sigma_values = np.asarray(sigma_values, dtype=float)
    total_sigma = np.sum(sigma_values)

    # 3. normalization
    #    avoid divide by zero
    if total_sigma > 1e-12:
      self.sigma_distribution = sigma_values / total_sigma
    else:
      # if sigma is zero everywhere
      self.sigma_distribution = np.zeros(self.n_nodes)

    # 4. Update Janus particle type
    if not np.allclose(self.sigma_distribution, self.sigma_distribution[0]):
      self.is_janus = True
    else:
      self.is_janus = False

  def get_surface_nodes(self, location = None, omega_axis_orientation = None, is_janus = False):
      """
      Calculates the positions of the surface nodes in the world frame.
      It rotates the body-frame nodes by the current orientation and then
      translates them to the body's current position.

      Returns:
          np.ndarray: An array of 3D vectors for each node's position in the world frame.
      """

      # Get location and orientation
      if location is None:
        location = self.location
      if omega_axis_orientation is None:
        omega_axis_orientation = self.omega_axis_orientation

      if is_janus:
        # Rotate each node from body frame to world frame using the quaternion
        rotation_matrix = omega_axis_orientation.rotation_matrix()
        rotated_nodes = np.dot(self.nodes_body_frame, rotation_matrix.T)
        # Translate rotated nodes to the body's position
        nodes_world = rotated_nodes + location
      else:
        # Translate original nodes to the body's position
        nodes_world = self.nodes_body_frame + location

      return nodes_world

  def update_v0_axis_from_omega_axis(self, omega_axis_orientation = None):
    """
    Initially, omega_axis corresponds to z-axis, v0_axis to x-axis.
    Re-calculates the orthogonal v_axis based on the current omega_axis.
    This is useful for initialization.
    """
    if omega_axis_orientation is None:
      omega_axis_orientation = self.omega_axis_orientation

    v0_axis_init = np.array([1.0, 0.0, 0.0])
    rotation_matrix = omega_axis_orientation.rotation_matrix()
    v0_axis = np.dot(v0_axis_init, rotation_matrix.T)

    self.v0_axis = v0_axis
    return v0_axis

  def update_omega_axis(self, omega_axis_orientation = None):
    """
    Initially, omega_axis corresponds to z-axis, v0_axis to x-axis.
    Re-calculates the orthogonal v_axis based on the current omega_axis.
    This is useful for initialization.
    """
    if omega_axis_orientation is None:
      omega_axis_orientation = self.omega_axis_orientation

    omega_axis_init = np.array([0.0, 0.0, 1.0])
    rotation_matrix = omega_axis_orientation.rotation_matrix()
    omega_axis = np.dot(omega_axis_init, rotation_matrix.T)

    self.omega_axis = omega_axis
    return omega_axis

  def calc_rot_matrix(self, location = None, orientation = None):
    '''
    Calculate the matrix R, where the i-th 3x3 block of R gives
    (R_i x) = -1 (r_i cross x).
    R has shape (3*N_nodes, 3).
    '''
    r_vectors = self.get_surface_nodes(location, orientation) - (self.location if location is None else location)
    rot_matrix = np.zeros((r_vectors.shape[0], 3, 3))
    rot_matrix[:, 0, 1] = r_vectors[:, 2]
    rot_matrix[:, 0, 2] = -r_vectors[:, 1]
    rot_matrix[:, 1, 0] = -r_vectors[:, 2]
    rot_matrix[:, 1, 2] = r_vectors[:, 0]
    rot_matrix[:, 2, 0] = r_vectors[:, 1]
    rot_matrix[:, 2, 1] = -r_vectors[:, 0]

    return np.reshape(rot_matrix, (3*self.n_nodes, 3))

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