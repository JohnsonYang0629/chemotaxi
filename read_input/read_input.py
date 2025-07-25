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

    self.mobility_alpha = float(self.options.get('mobility_alpha') or 1)
    self.radius = float(self.options.get('radius') or 1.0)
    self.intrinsic_linear_velocity = float(self.options.get('intrinsic_linear_velocity') or 1.0)
    self.intrinsic_angular_velocity = float(self.options.get('intrinsic_angular_velocity') or 1.0)
    self.emission_rate = float(self.options.get('emission_rate') or 1.0)
    self.persistence_length = float(self.options.get('persistence_length') or 1.0)
    self.gamma_t = float(self.options.get('translational_noise_gamma') or 500.0)
    self.gamma_r = float(self.options.get('rotational_noise_gamma') or 500.0)
    self.peclet_number = float(self.options.get('peclet_number') or 1.0)

    self.initial_position_2D = np.fromstring(self.options.get('initial_position_2D') or '0 0', sep=' ')
    self.initial_orientation_2D_vector = np.fromstring(self.options.get('initial_orientation_2D_vector') or '0 0',
                                                       sep=' ')

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

    self.initial_position_3D = np.fromstring(self.options.get('initial_position_3D') or '0 0 0', sep=' ')
    # Prepare quaternion for omega axis
    orientation_quaternion_input = np.fromstring(self.options.get('initial_orientation_3D_quaternion') or
                                                 '1 0 0 0', sep=' ')
    orientation = [float(orientation_quaternion_input[0]), float(orientation_quaternion_input[1]),
                   float(orientation_quaternion_input[2]), float(orientation_quaternion_input[3])]
    norm_orientation = np.linalg.norm(orientation)
    orientation_quaternion = Quaternion(orientation / norm_orientation)
    self.initial_orientation_3D_quaternion = orientation_quaternion

    return
