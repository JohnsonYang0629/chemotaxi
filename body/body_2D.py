import numpy as np
import copy
import sys


class Body2D(object):
  '''
  Small class to handle a single body in 2D domain.
  '''
  def __init__(self, location, orientation, structure_ref_config, n_steps):
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
    self.orientation_history = np.zeros([n_steps + 1, 2])
    # Reference configuration. Coordinates of droplet for polar coordinate theta [0]
    # and location = np.array[0, 0]) as a np.array.shape = (1, 2)
    # self.rotation_matrix = None
    # Some default functions
    self.function_force = self.default_none
    self.function_torque = self.default_none
    self.prescribed_velocity = np.array([0.0, 0.0, 0.0])
    self.velocities_previous_step = np.array([0.0, 0.0, 0.0])
    self.chem_surface_gradient = np.array([0.0, 0.0])

    # self.structure_ref_config = np.copy(structure_ref_config)
    self.nodes_body_frame = np.copy(structure_ref_config)
    self.n_nodes = len(self.nodes_body_frame)

    self.sigma_distribution = np.ones(self.n_nodes) / self.n_nodes
    self.is_janus = False
    self.chem_surface_gradient = np.zeros(2)
    self.chem_torque_gradient = 0.0

    self.peclet_number = 1.0
    self.mobility_alpha = 1.0
    self.mobility_distribution = np.ones(self.n_nodes) * self.mobility_alpha
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

  def set_sigma_distribution(self, sigma_values):
    """
    Set the sigma (chemical release/absorption rate) distribution for 2D particles.
    Args:
        sigma_values (list or np.ndarray): Each node's corresponding sigma value
                                           (e.g., -1.0 for sink, 1.0 for source).
    """
    # 1. Check input validity to ensure the array matches the mesh size
    if sigma_values is None or len(sigma_values) != self.n_nodes:
      print("Warning: Input sigma_values invalid, default uniform distribution will be used.")
      self.sigma_distribution = np.ones(self.n_nodes) / self.n_nodes
      self.is_janus = False
      return

    # 2. Read the user-defined surface release rates and convert to a flat numpy array
    sigma_values = np.asarray(sigma_values, dtype=float).flatten()

    # 3. Scale by the number of nodes (1/N) instead of normalizing by the sum.
    # This preserves the net physical release/absorption rates and avoids division
    # by zero when sources and sinks perfectly cancel each other out (sum = 0).
    self.sigma_distribution = sigma_values / self.n_nodes

    # 4. Determine if it is a Janus particle.
    # If the surface distribution is not perfectly uniform, treat it as a Janus particle.
    if not np.allclose(self.sigma_distribution, self.sigma_distribution[0]):
      self.is_janus = True
    else:
      self.is_janus = False

  def set_sigma_distribution_test(self, sigma_values):
      """
      Set and normalize sigma distribution for 2D particles.
      Args:
          sigma_values (list or np.ndarray): each node's corresponding sigma value.
      """
      # 1. check the length of sigma value is identical to the nodes.
      if sigma_values is None or len(sigma_values) != self.n_nodes:
        print("Warning: Input sigma_values invalid, default uniform distribution will be used.")
        # if costume distribution invalid, we use non-Janus particle.
        self.sigma_distribution = np.ones(self.n_nodes) / self.n_nodes
        self.is_janus = False
        return

      # 2. normalize with the sum of total sigma
      sigma_values = np.asarray(sigma_values, dtype=float).flatten()
      total_sigma = np.sum(sigma_values)

      # 3. normalization avoid divide by zero
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

  def set_mobility_distribution(self, mobility_values):
    """
    Set the local mobility (alpha) distribution for 2D particles.
    Args:
        mobility_values (list or np.ndarray): Each node's corresponding mobility value.
    """
    # 1. Check validity. If invalid or not provided, keep the uniform global alpha.
    if mobility_values is None or len(mobility_values) != self.n_nodes:
      print(
        f"Info: No valid mobility distribution provided for Particle {self.ID}. Using uniform scalar alpha = {self.mobility_alpha}.")
      self.mobility_distribution = np.ones(self.n_nodes) * self.mobility_alpha
      return
    # 2. Read the user-defined surface mobility rates and convert to a flat numpy array
    self.mobility_distribution = np.asarray(mobility_values, dtype=float).flatten()

  def get_surface_nodes(self, location=None, orientation=None, is_janus=False):
    """
    Calculates the positions of the 2D surface nodes in the world (lab) frame.
    orientation expects a 2D unit vector [x, y].
    """
    if location is None:
      location = self.location
    if orientation is None:
      orientation = self.orientation

    # In 2D, the unit vector orientation represents [cos(theta), sin(theta)]
    # As long as it's a Janus particle or has rotated from [1.0, 0.0], perform rotation
    if is_janus or not (np.isclose(orientation[0], 1.0) and np.isclose(orientation[1], 0.0)):
      cos_t = orientation[0]
      sin_t = orientation[1]

      # 2D rotation matrix
      rot_matrix = np.array([[cos_t, -sin_t],
                             [sin_t, cos_t]])

      # Translate the nodes after rotating them from the body frame
      rotated_nodes = np.dot(self.nodes_body_frame, rot_matrix.T)
      nodes_world = rotated_nodes + location
    else:
      # Direct translation without rotation
      nodes_world = self.nodes_body_frame + location

    return nodes_world

  def default_none(self, *args, **kwargs):
    return None
