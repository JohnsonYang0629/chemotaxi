import numpy as np
from functools import partial


class ChemoIntegrator2D(object):

    def __init__(self, bodies, scheme, domain, numerical_method):
        """
        Initialize object
        """
        self.bodies = bodies
        self.scheme = scheme
        self.domain = domain
        self.numerical_method = numerical_method

        # Other variables
        self.velocities = None
        self.velocities_previous_step = None
        self.first_step = True
        self.intrinsic_velocity = np.array([0, 0])
        self.gamma_r = 0.0
        self.gamma_t = 0.0

        # Optional variables
        self.calc_surface_gradient_circle_numba_optimized = None
        self.rotation_matrix_2d = None

    def advance_time_step(self, dt, *args, **kwargs):
        """
        Main entry point for advancing time.
        This method now iterates through all the bodies it owns and updates them one by one.
        """
        # A deep copy of states might be needed for more complex integration schemes
        # to prevent using partially updated states within the same time step.
        # For this forward scheme, iterating is sufficient.
        for body_to_update in self.bodies:
            # Dispatch to the appropriate scheme for the current body.
            getattr(self, self.scheme)(body_to_update, dt, *args, **kwargs)

    def history_local_compose_2d(self, body_to_update, dt, *args, **kwargs):
        """
        History part:
        [0, (N-1)dt]
        using the trapezoidal integration scheme.
        Local part:
        [(N-1)dt, Ndt]
        """
        while True:
            step = kwargs.get('step')
            peclet_number = body_to_update.peclet_number
            mobility_alpha = body_to_update.mobility_alpha

            if self.first_step == False:
                # Use history-local compose method
                chem_force = self.calc_surface_gradient_circle_numba_optimized(body_to_update, self.bodies, *args, **kwargs)
                chem_prop = mobility_alpha/(2 * np.pi) * chem_force
                angular_velocity = self.intrinsic_velocity[1]
                # Two-step Adams-Bashforth method
                if self.numerical_method == "adams_bashforth_2":
                    angular_velocity_dt = (1.5 * angular_velocity - 0.5 * body_to_update.velocities_previous_step[2]) * dt
                # Forward Euler method
                if self.numerical_method == "forward_euler":
                    angular_velocity_dt = angular_velocity * dt
                # Stochastic First Order
                if self.numerical_method == "stochastic_first_order":
                    random_rotation = np.random.randn()
                    rotational_noise_term = np.sqrt(2 / self.gamma_r) * random_rotation * np.sqrt(dt)
                    angular_velocity_dt = angular_velocity * dt + rotational_noise_term

                orientation_new = np.dot(self.rotation_matrix_2d(angular_velocity_dt), body_to_update.orientation)
                body_to_update.orientation = orientation_new
                linear_velocity_compose = self.intrinsic_velocity[0] * orientation_new + chem_prop
                # Two-step Adams-Bashforth method
                if self.numerical_method == "adams_bashforth_2":
                    location_new = body_to_update.location + (1.5 * linear_velocity_compose - 0.5 * body_to_update.velocities_previous_step[0:2]) * dt
                # Forward Euler method
                if self.numerical_method == "forward_euler":
                    location_new = body_to_update.location + linear_velocity_compose * dt
                # Stochastic First Order
                if self.numerical_method == "stochastic_first_order":
                    random_translation = np.random.randn(2)
                    translational_noise_term = np.sqrt(2 / self.gamma_t) * random_translation * np.sqrt(dt)
                    location_new = body_to_update.location + linear_velocity_compose * dt + translational_noise_term

                body_to_update.location = location_new
                velocity = np.append(linear_velocity_compose, angular_velocity)
                body_to_update.prescribed_velocity = velocity
                body_to_update.velocities_previous_step = velocity
                body_to_update.chem_surface_gradient = chem_force

            else:
                # Use forward Euler method for the first step
                chem_force = self.calc_surface_gradient_circle_numba_optimized(body_to_update, self.bodies, *args, **kwargs)
                chem_prop = mobility_alpha/(2 * np.pi) * chem_force
                angular_velocity = self.intrinsic_velocity[1]  # noise required
                angular_velocity_dt = angular_velocity * dt
                orientation_new = np.dot(self.rotation_matrix_2d(angular_velocity_dt), body_to_update.orientation)
                body_to_update.orientation = orientation_new
                linear_velocity_compose = self.intrinsic_velocity[0] * orientation_new + chem_prop
                location_new = body_to_update.location + linear_velocity_compose[0:2] * dt
                body_to_update.location = location_new
                velocity = np.append(linear_velocity_compose, angular_velocity)
                body_to_update.prescribed_velocity = velocity
                body_to_update.chem_surface_gradient = chem_force

            # Update configuration
            body_to_update.location_history[step + 1, :] = location_new
            self.first_step = False
            body_to_update.velocities_previous_step = velocity

            return
