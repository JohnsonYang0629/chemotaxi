import numpy as np
from functools import partial


class ChemoIntegrator2D(object):

    def __init__(self, body, scheme, domain, numerical_method):
        """
        Initialize object
        """
        self.body = body
        self.scheme = scheme
        self.domain = domain
        self.numerical_method = numerical_method

        # Other variables
        self.velocities = None
        self.velocities_previous_step = None
        self.first_step = True
        self.peclet_number = 0.0
        self.mobility_alpha = 0.0
        self.intrinsic_velocity = np.array([0, 0])

        # Optional variables
        self.calc_surface_gradient_circle = None
        self.rotation_matrix_2d = None

    def advance_time_step(self, dt, *args, **kwargs):
        """
        Advance time step with key: integrator self.scheme
        """
        return getattr(self, self.scheme)(dt, *args, **kwargs)

    def history_local_compose_2d(self, dt, *args, **kwargs):
        """
        History part:
        [0, (N-1)dt]
        using the trapezoidal integration scheme.
        Local part:
        [(N-1)dt, Ndt]
        """
        while True:
            step = kwargs.get('step')
            body = self.body
            if self.first_step == False:
                # Use history-local compose method
                chem_force, chem_local, chem_history = self.calc_surface_gradient_circle(self.body, *args, **kwargs)
                chem_prop = self.mobility_alpha/(2 * np.pi) * chem_force
                angular_velocity = self.intrinsic_velocity[1]  # noise required
                # Two-step Adams-Bashforth method
                if self.numerical_method == "adams_bashforth_2":
                    angular_velocity_dt = (1.5 * angular_velocity - 0.5 * self.velocities_previous_step[2]) * dt
                # Forward Euler method
                if self.numerical_method == "forward_euler":
                    angular_velocity_dt = angular_velocity * dt
                orientation_new = np.dot(self.rotation_matrix_2d(angular_velocity_dt), body.orientation)
                body.orientation = orientation_new
                linear_velocity_compose = orientation_new + chem_prop
                # Two-step Adams-Bashforth method
                if self.numerical_method == "adams_bashforth_2":
                    location_new = body.location + (1.5 * linear_velocity_compose - 0.5 * self.velocities_previous_step[0:2]) * dt
                # Forward Euler method
                if self.numerical_method == "forward_euler":
                    location_new = body.location + linear_velocity_compose * dt
                body.location = location_new
                velocity = np.append(linear_velocity_compose, angular_velocity)
                body.prescribed_velocity = velocity
                body.chem_surface_gradient = chem_force
                body.grad_chem_local_part = chem_local
                body.grad_chem_history_part = chem_history

            else:
                # Use forward Euler method for the first step
                chem_force, chem_local, chem_history = self.calc_surface_gradient_circle(self.body, *args, **kwargs)
                chem_prop = self.mobility_alpha/(2 * np.pi) * chem_force
                angular_velocity = self.intrinsic_velocity[1]  # noise required
                angular_velocity_dt = angular_velocity * dt
                orientation_new = np.dot(self.rotation_matrix_2d(angular_velocity_dt), body.orientation)
                body.orientation = orientation_new
                linear_velocity_compose = orientation_new + chem_prop
                location_new = body.location + linear_velocity_compose[0:2] * dt
                body.location = location_new
                velocity = np.append(linear_velocity_compose, angular_velocity)
                body.prescribed_velocity = velocity
                body.chem_surface_gradient = chem_force
                body.grad_chem_local_part = chem_local
                body.grad_chem_history_part = chem_history

            # Update configuration
            body.location_history[step + 1, :] = location_new
            self.first_step = False
            self.velocities_previous_step = velocity


            return
