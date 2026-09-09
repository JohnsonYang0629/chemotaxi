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
        self.history_local_compose_2d_multi_body = None
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
        Updates the state of a single 2D body.
        """
        step = kwargs.get('step')
        force_grad = np.zeros(2)
        torque_grad = 0.0

        if not self.first_step:
            force_grad, torque_grad = self.history_local_compose_2d_multi_body(
                target_body=body_to_update,
                all_bodies=self.bodies,
                dt=dt,
                step=step
            )

        # chem_prop = -(body_to_update.mobility_alpha / (2 * np.pi)) * force_grad
        # chem_prop = -(1.0 / (2 * np.pi)) * force_grad
        chem_prop = -1.0 * force_grad

        if body_to_update.is_janus:
            # Janus particle: turn off omega_0, compute torque induced by asymmetric chemical field
            # mobility_rotational = body_to_update.mobility_alpha * 0.75
            # chem_torque = -(mobility_rotational / (2 * np.pi)) * torque_grad
            # chem_torque = -(0.75 / (2 * np.pi)) * torque_grad
            chem_torque = -0.75 * torque_grad
            angular_velocity = self.intrinsic_velocity[1] + chem_torque
            # In 2D, the orientation unit vector is exactly the v0_axis
            intrinsic_swim_velocity = body_to_update.orientation * self.intrinsic_velocity[0]
        else:
            # Isotropic particle: keep intrinsic chirality omega_0
            chem_torque = 0.0
            angular_velocity = self.intrinsic_velocity[1] + chem_torque
            intrinsic_swim_velocity = body_to_update.orientation * self.intrinsic_velocity[0]

        linear_velocity_compose = intrinsic_swim_velocity + chem_prop

        # --- Update vectors using 2D rotation matrix ---
        if self.numerical_method == "adams_bashforth_2" and not self.first_step:
            # Position update
            location_new = body_to_update.location + (
                        1.5 * linear_velocity_compose - 0.5 * body_to_update.velocities_previous_step[0:2]) * dt

            # Angular velocity AB2 prediction
            angular_velocity_ab2 = 1.5 * angular_velocity - 0.5 * body_to_update.velocities_previous_step[2]

            # Vector rotation
            rot_matrix = self.rotation_matrix_2d(angular_velocity_ab2 * dt)
            orientation_new = np.dot(rot_matrix, body_to_update.orientation)
        else:
            # Default Forward Euler
            location_new = body_to_update.location + linear_velocity_compose * dt

            rot_matrix = self.rotation_matrix_2d(angular_velocity * dt)
            orientation_new = np.dot(rot_matrix, body_to_update.orientation)

        # Normalize to prevent vector length drift due to numerical rounding
        orientation_new = orientation_new / np.linalg.norm(orientation_new)

        # Update states
        body_to_update.location = location_new
        body_to_update.orientation = orientation_new

        body_to_update.location_history[step + 1, :] = location_new
        body_to_update.orientation_history[step + 1, :] = orientation_new

        # Record prescribed velocity, in 2D it is [v_x, v_y, omega]
        body_to_update.prescribed_velocity = np.array(
            [linear_velocity_compose[0], linear_velocity_compose[1], angular_velocity])
        body_to_update.chem_surface_gradient = force_grad
        body_to_update.chem_torque_gradient = torque_grad
        body_to_update.velocities_previous_step = body_to_update.prescribed_velocity

        self.first_step = False

    def history_local_compose_2d_nonjanus(self, body_to_update, dt, *args, **kwargs):
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
