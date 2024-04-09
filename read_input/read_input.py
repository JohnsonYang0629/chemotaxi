'''
Simple class to read the input files to run a simulation.
'''

import numpy as np
import ntpath
import sys


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
    self.scheme = str(self.options.get('scheme') or 'HLP_2D')

    self.mobility_alpha = int(self.options.get('mobility_alpha') or 1)
    self.radius = float(self.options.get('radius') or 1.0)
    self.intrinsic_linear_velocity = float(self.options.get('intrinsic_linear_velocity') or 1.0)
    self.intrinsic_angular_velocity = float(self.options.get('intrinsic_angular_velocity') or 1.0)
    self.emission_rate = float(self.options.get('emission_rate') or 1.0)
    self.persistence_length = float(self.options.get('persistence_length') or 1.0)
    self.peclet_number = float(self.options.get('peclet_number') or 1.0)

    self.initial_position = np.fromstring(self.options.get('initial_position') or '0 0', sep=' ')
    self.initial_orientation = np.fromstring(self.options.get('initial_orientation') or '0 0 0 0', sep=' ')

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

    return
