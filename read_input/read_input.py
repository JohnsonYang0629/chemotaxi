'''
Simple class to read the input files to run a simulation.
'''

import numpy as np
import ntpath
import sys
from integrator.quaternion import Quaternion


class ReadInput(object):
  """
  Simple class to read the input files to generate server required input/command files.
  """

  def __init__(self, entries):
    """ Constructor takes the name of the input file """
    self.entries = entries
    self.input_file = entries
    self.options = {}
    number_of_structures = 0
    number_of_chem_dists = 0
    number_of_mob_dists = 0

    # Read input file
    comment_symbols = ['#']   
    with open(self.input_file, 'r') as f:
      # Loop over lines
      for line in f:
        # Strip comments
        if comment_symbols[0] in line:
          line, comment = line.split(comment_symbols[0], 1)

        # Save options to dictionary, Value may be more than one word
        line = line.strip()
        if line != '':
          option, value = line.split(None, 1)
          if option == 'structure':
            option += str(number_of_structures)
            number_of_structures += 1
          elif option == 'chemical_distribution':
            option += str(number_of_chem_dists)
            number_of_chem_dists += 1
          elif option == 'mobility_distribution':
            option += str(number_of_mob_dists)
            number_of_mob_dists += 1
          self.options[option] = value

    # Set options to test or default values
    self.job_name = str(self.options.get('job_name') or 'None')
    self.job_type = str(self.options.get('job_type') or 'None')
    self.domain = str(self.options.get('domain') or '2D')
    self.scheme = str(self.options.get('scheme') or 'history_local_compose_2d')
    self.particle_type = str(self.options.get('particle_type') or 'non_janus')
    self.acceleration = str(self.options.get('acceleration') or 'numba')
    self.core = int(self.options.get('core') or 1)
    self.numerical_method = str(self.options.get('numerical_method') or 'forward_euler')

    self.droplet_num = int(self.options.get('droplet_num') or 1)

    self.radius = float(self.options.get('radius') or 1.0)
    self.intrinsic_linear_velocity = float(self.options.get('intrinsic_linear_velocity') or 1.0)
    self.intrinsic_angular_velocity = float(self.options.get('intrinsic_angular_velocity') or 1.0)
    self.emission_rate = float(self.options.get('emission_rate') or 1.0)
    self.persistence_length = float(self.options.get('persistence_length') or 1.0)
    self.gamma_t = float(self.options.get('translational_noise_gamma') or 500.0)
    self.gamma_r = float(self.options.get('rotational_noise_gamma') or 500.0)

    self.surface_disc_num = int(self.options.get('surface_disc_num') or 2)

    self.dt = float(self.options.get('dt') or 0.0)
    self.n_steps = int(self.options.get('n_steps') or 1)
    self.n_save = int(self.options.get('n_save') or 1)
    self.initial_step = int(self.options.get('initial_step') or 0)

    self.kT = float(self.options.get('kT') or 1.0)
    self.eta = float(self.options.get('eta') or 1.0)
    self.g = float(self.options.get('g') or 1.0)

    self.output_name = str(self.options.get('output_name') or 'run')
    self.save_clones = str(self.options.get('save_clones') or 'one_file')
    self.structure = str.split(str(self.options.get('structure0')))
    self.chemical_distribution_file = str.split(str(self.options.get('chemical_distribution') or 'None'))
    self.random_state = str(self.options.get('random_state') or 'None')

    self.mobility_alphas = np.fromstring(self.options.get('mobility_alpha') or '1.0', sep=' ')
    self.peclet_numbers = np.fromstring(self.options.get('peclet_number') or '1.0', sep=' ')

    self.particle_types = str.split(self.options.get('particle_type') or 'non_janus')

    if len(self.particle_types) != self.droplet_num:
      print(f"Warning: Number of particle_type values does not match droplet_num. Using first value for all.")
      self.particle_types = [self.particle_types[0]] * self.droplet_num
    if len(self.mobility_alphas) != self.droplet_num:
      print(f"Warning: Number of mobility_alpha values does not match droplet_num. Using first value for all.")
      self.mobility_alphas = np.full(self.droplet_num, self.mobility_alphas[0])
    if len(self.peclet_numbers) != self.droplet_num:
      print(f"Warning: Number of peclet_number values does not match droplet_num. Using first value for all.")
      self.peclet_numbers = np.full(self.droplet_num, self.peclet_numbers[0])

    if self.domain == "2D":
      # For 2D
      # Read the long string of 2D positions and reshape it into a (droplet_num, 2) array.
      initial_pos_2d_flat = np.fromstring(self.options.get('initial_position_2D') or '0 0', sep=' ')
      if len(initial_pos_2d_flat) != self.droplet_num * 2:
        sys.exit(
          f"Error: Mismatch between droplet_num ({self.droplet_num}) and number of initial_position_2D coordinates.")
      self.initial_positions_2D = initial_pos_2d_flat.reshape((self.droplet_num, 2))

      # Read the long string of 2D orientations and reshape it.
      initial_orient_2d_flat = np.fromstring(self.options.get('initial_orientation_2D_vector') or '1 0', sep=' ')
      if len(initial_orient_2d_flat) != self.droplet_num * 2:
        sys.exit(
          f"Error: Mismatch between droplet_num ({self.droplet_num}) and number of initial_orientation_2D_vector coordinates.")

      # Ensure the input 2D orientation is normalized to a unit vector [x, y]
      orient_vectors_2d = initial_orient_2d_flat.reshape((self.droplet_num, 2))
      norms = np.linalg.norm(orient_vectors_2d, axis=1, keepdims=True)
      # Avoid division by zero for zero vectors
      self.initial_orientations_2D = np.where(norms > 1e-9, orient_vectors_2d / norms, orient_vectors_2d)

    elif self.domain == "3D":
      # For 3D
      # Read the long string of 3D positions and reshape it into a (droplet_num, 3) array.
      initial_pos_3d_flat = np.fromstring(self.options.get('initial_position_3D') or '0 0 0', sep=' ')
      if len(initial_pos_3d_flat) != self.droplet_num * 3:
        sys.exit(
          f"Error: Mismatch between droplet_num ({self.droplet_num}) and number of initial_position_3D coordinates.")
      self.initial_positions_3D = initial_pos_3d_flat.reshape((self.droplet_num, 3))

      # Read the long string of 3D quaternion orientations and reshape it.
      initial_orient_3d_flat = np.fromstring(self.options.get('initial_orientation_3D_quaternion') or '1 0 0 0',
                                             sep=' ')
      if len(initial_orient_3d_flat) != self.droplet_num * 4:
        sys.exit(
          f"Error: Mismatch between droplet_num ({self.droplet_num}) and number of initial_orientation_3D_quaternion coordinates.")
      orientations_3d_reshaped = initial_orient_3d_flat.reshape((self.droplet_num, 4))

      # Create a list of Quaternion objects.
      self.initial_orientations_3D_quaternion = []
      for orient_array in orientations_3d_reshaped:
        norm = np.linalg.norm(orient_array)
        # Normalize the quaternion to ensure it represents a pure rotation.
        q_values = orient_array / norm if norm > 1e-9 else orient_array
        self.initial_orientations_3D_quaternion.append(Quaternion(q_values))

    self.structures = []
    for i in range(number_of_structures):
      self.structures.append(str.split(str(self.options.get(f'structure{i}'))))

    self.chemical_distribution_files = []
    for i in range(number_of_chem_dists):
      self.chemical_distribution_files.append(str.split(str(self.options.get(f'chemical_distribution{i}')))[0])

    self.mobility_distribution_files = []
    for i in range(number_of_mob_dists):
      self.mobility_distribution_files.append(str.split(str(self.options.get(f'mobility_distribution{i}')))[0])

    return
