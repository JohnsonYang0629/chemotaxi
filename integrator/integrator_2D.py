import numpy as np


class ChemoIntegrator2D(object):

  def __init__(self, body, scheme, domain):
    '''
    Initialize object
    '''
    self.body = body
    self.scheme = scheme

    # Other variables
    self.domain = domain
    self.eta = None
    self.a = None
    self.velocities = None
    self.velocities_previous_step = None
    self.first_step = True
    self.kT = 0.0
    self.tolerance = 1e-08
    self.rf_delta = 1e-03
    self.invalid_configuration_count = 0
    self.det_iterations_count = 0
    self.stoch_iterations_count = 0


    # Optional variables
    self.periodic_length = None
    self.calc_slip = None
    self.calc_force_torque = None
    self.mobility_inv_blobs = None
    self.first_guess = None
    self.mobility_vector_prod = None

  def advance_time_step(self, dt, *args, **kwargs):
    '''
    Advance time step with integrator self.scheme
    '''
    return getattr(self, self.scheme)(dt, *args, **kwargs)

  def history_local_compose_2d(self, dt, *args, **kwargs):
    '''
    Take a time step of length dt using the trapezoidal integration scheme.
    '''
    while True:
      # Call preprocess
      preprocess_result = self.preprocess(self.bodies)

      fix_pos = self.fix_pos
      fix_height = self.fix_height

      # Solve mobility problem
      sol_precond = self.solve_mobility_problem(x0 = self.first_guess, save_first_guess = True, update_PC = self.update_PC, step = kwargs.get('step'))

      # Extract velocities
      velocities = np.reshape(sol_precond[3*self.Nblobs: 3*self.Nblobs + 6*len(self.bodies)], (len(self.bodies) * 6))
      # Update location and orientation
      if self.first_step == False:
        # Use Adams-Bashforth
        for k, b in enumerate(self.bodies):
          b.location_new = b.location + (1.5 * velocities[6*k:6*k+3] - 0.5 * self.velocities_previous_step[6*k:6*k+3]) * dt
          if fix_pos == 'True':
            b.location_new = b.location   # Fixed position
            print("fixed location = ", b.location_new)
          elif fix_height == 'True':
            b.location_new[2] = b.location[2]
            print("fixed height = ", b.location_new)
          else:
            print("moving location = ", b.location_new)
          quaternion_dt = Quaternion.from_rotation((1.5 * velocities[6*k+3:6*k+6] - 0.5 * self.velocities_previous_step[6*k+3:6*k+6]) * dt)
          b.orientation_new = quaternion_dt * b.orientation
          b.prescribed_velocity = velocities[6*k:6*k+6]

          # print("orientation = ", b.orientation_new)
          # print("velocity =", b.prescribed_velocity)
      else:
        # Use forward Euler
        for k, b in enumerate(self.bodies):
          b.location_new = b.location + velocities[6*k:6*k+3] * dt
          if fix_pos == 'True':
            b.location_new = b.location   # Fixed position
            print("fixed location = ", b.location_new)
          elif fix_height == 'True':
            b.location_new[2] = b.location[2]
            print("fixed height = ", b.location_new)
          else:
            print("moving location = ", b.location_new)
          quaternion_dt = Quaternion.from_rotation((velocities[6*k+3:6*k+6]) * dt)
          b.orientation_new = quaternion_dt * b.orientation
          b.prescribed_velocity = velocities[6*k:6*k+6]

          # print("orientation = ", b.orientation_new)
          # print("velocity =", b.prescribed_velocity)

      # Call postprocess
      postprocess_result = self.postprocess(self.bodies)

      # Check positions, if valid return
      if self.check_positions(new = 'new', old = 'current', update_in_success = True, domain = self.domain) is True:
        # Save velocities for next step
        self.first_step = False
        self.velocities_previous_step = velocities
        return

    return