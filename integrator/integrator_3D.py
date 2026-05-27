import numpy as np
from integrator.quaternion import Quaternion


class ChemoIntegrator3D(object):

    def __init__(self, bodies, scheme, domain, numerical_method):
        """
        Initialize object
        """
        self.bodies = bodies
        self.scheme = scheme
        self.domain = domain
        self.numerical_method = numerical_method

        # Other variables
        self.velocities = None  # [vx, vy, vz, wx, wy, wz]
        self.velocities_previous_step = None
        self.first_step = True
        self.intrinsic_velocity = np.array([0, 0])  # compact vector [v_0,omega_0]
        self.gamma_r = 0.0  # Rotational diffusion
        self.gamma_t = 0.0  # Translational diffusion

        # Optional variables
        self.calc_tangential_grad_3D = None
        self.history_local_compose_3d_distribution = None
        self.history_local_compose_3d_point = None
        self.history_local_compose_3d_multi_body = None
        self.rotation_matrix_3d = None

    def advance_time_step(self, dt, *args, **kwargs):
        """
        Main entry point for advancing the simulation by one time step.
        This method iterates through all bodies and calls the designated scheme
        function to update each one individually.
        """
        # The scheme (e.g., 'history_local_compose_3d') is called for each body.
        for body_to_update in self.bodies:
            getattr(self, self.scheme)(body_to_update, dt, *args, **kwargs)

    def noise(self, dt, *args, **kwargs):
        return

    def history_local_compose_3d(self, body_to_update, dt, *args, **kwargs):
        """
        Updates the state of a single 3D body (`body_to_update`) for one time step.
        It computes the chemical interactions from ALL bodies in the simulation
        to determine the macroscopic chemotactic force and torque on this specific body.
        """
        while True:
            step = kwargs.get('step')
            force_grad = np.zeros(3)
            torque_grad = np.zeros(3)

            # On the first step, we assume zero initial chemical gradient.
            # For subsequent steps, we calculate the full multi-body interactions.
            if not self.first_step:
                force_grad, torque_grad = self.history_local_compose_3d_multi_body(
                    target_body=body_to_update,
                    all_bodies=self.bodies,
                    dt=dt,
                    step=step
                )

            # Calculate translational chemical propulsion (Fc)
            chem_prop = -(body_to_update.mobility_alpha / (4 * np.pi)) * force_grad

            # Determine rotational dynamics based on particle type
            if body_to_update.is_janus:
                # Janus particle: No intrinsic structural rotation (omega_0 = 0)
                # Rotational mobility (\Lambda_r) is 3/4 of translational mobility (\Lambda) for a solid sphere
                mobility_rotational = body_to_update.mobility_alpha * 0.75
                chem_torque = -(mobility_rotational / (4 * np.pi)) * torque_grad

                angular_velocity_vector = chem_torque
                intrinsic_swim_velocity = body_to_update.v0_axis * self.intrinsic_velocity[0]
            else:
                # Isotropic droplet: Driven by intrinsic chirality, auto-chemotactic macroscopic torque is negligible
                chem_torque = np.zeros(3)
                angular_velocity_vector = self.intrinsic_velocity[1] * body_to_update.omega_axis + chem_torque
                intrinsic_swim_velocity = body_to_update.v0_axis * self.intrinsic_velocity[0]

            # Compose the final linear velocity
            linear_velocity_compose = intrinsic_swim_velocity + chem_prop

            # --- Update position and orientation based on the chosen numerical method ---

            if self.numerical_method == "adams_bashforth_2" and not self.first_step:
                # Two-step Adams-Bashforth method (AB2)
                # velocities_previous_step stores [vx, vy, vz, wx, wy, wz]
                location_new = body_to_update.location + (
                            1.5 * linear_velocity_compose - 0.5 * body_to_update.velocities_previous_step[0:3]) * dt

                omega_axis_quaternion_dt = Quaternion.from_rotation(
                    (1.5 * angular_velocity_vector - 0.5 * body_to_update.velocities_previous_step[3:6]) * dt)
                body_to_update.omega_axis_orientation = omega_axis_quaternion_dt * body_to_update.omega_axis_orientation
                body_to_update.location = location_new

            elif self.numerical_method == "stochastic_first_order":
                # Add stochastic terms for Brownian motion
                random_rotation_vec = np.random.randn(3)
                stochastic_rotation_term = np.sqrt(2 / self.gamma_r) * random_rotation_vec * np.sqrt(dt)

                random_translation_vec = np.random.randn(3)
                stochastic_translation_term = np.sqrt(2 / self.gamma_t) * random_translation_vec * np.sqrt(dt)

                # Update orientation with deterministic and stochastic parts
                rotation_increment = Quaternion.from_rotation(angular_velocity_vector * dt + stochastic_rotation_term)
                body_to_update.omega_axis_orientation = rotation_increment * body_to_update.omega_axis_orientation

                # Update translational position
                body_to_update.location += linear_velocity_compose * dt + stochastic_translation_term

            else:
                # Default to Forward Euler method (also used for the first step of AB2 initialization)
                rotation_increment = Quaternion.from_rotation(angular_velocity_vector * dt)
                body_to_update.omega_axis_orientation = rotation_increment * body_to_update.omega_axis_orientation
                body_to_update.location += linear_velocity_compose * dt

            # --- Finalize state update for the current body ---

            # Update directional axes based on the new orientation
            body_to_update.update_omega_axis()
            body_to_update.update_v0_axis_from_omega_axis()

            # Store history for the next time step's calculations
            body_to_update.location_history[step + 1, :] = body_to_update.location
            # Save orientation in [x, y, z, w] format
            body_to_update.orientation = body_to_update.omega_axis_orientation.flip_self()
            body_to_update.orientation_history[step + 1, :] = body_to_update.orientation

            # Store velocity and gradient for output and potential use in higher-order integrators
            body_to_update.prescribed_velocity = np.append(linear_velocity_compose, angular_velocity_vector)
            body_to_update.chem_surface_gradient = force_grad
            body_to_update.chem_torque_gradient = torque_grad
            body_to_update.velocities_previous_step = body_to_update.prescribed_velocity

            # Flag off first step after initialization
            self.first_step = False
            return

    def history_local_compose_3d_single_body(self, dt, *args, **kwargs):
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
            force_grad = np.zeros(3)
            torque_grad = np.zeros(3)  # Torque is zero by default
            chem_torque = np.zeros(3)  # Torque is zero by default

            if self.first_step == False:
                # Use history-local compose method
                # chem_force = self.calc_tangential_grad_3D(self.body, *args, **kwargs)
                # chem_prop = self.mobility_alpha / (4 * np.pi) * chem_force

                if body.is_janus:
                    # Distribution case: returns force and torque gradients
                    force_grad, torque_grad = self.history_local_compose_3d_distribution(self.body, *args, **kwargs)
                    chem_prop = (self.mobility_alpha / (4 * np.pi)) * force_grad
                    chem_torque = (self.mobility_alpha / (4 * np.pi)) * torque_grad
                else:
                    # Point case: only returns a force gradient
                    force_grad = self.history_local_compose_3d_point(self.body, *args, **kwargs)
                    chem_prop = (self.mobility_alpha / (4 * np.pi)) * force_grad

                intrinsic_swim_velocity = body.v0_axis * self.intrinsic_velocity[0]
                linear_velocity_compose = intrinsic_swim_velocity + chem_prop

                omega_axis = body.update_omega_axis(body.omega_axis_orientation)
                angular_velocity_vector = self.intrinsic_velocity[1] * omega_axis + chem_torque

                # Two-step Adams-Bashforth method
                if self.numerical_method == "adams_bashforth_2":
                    location_new = body.location \
                                   + (1.5 * linear_velocity_compose - 0.5 * self.velocities_previous_step[0:3]) * dt

                    omega_axis_quaternion_dt = Quaternion.from_rotation((1.5 * angular_velocity_vector -
                                                                         0.5 * self.velocities_previous_step[4:6]) * dt)
                    body.location = location_new
                    body.omega_axis_orientation = omega_axis_quaternion_dt * body.omega_axis_orientation
                    body.v0_axis = body.update_v0_axis_from_omega_axis(body.omega_axis_orientation)
                    velocity = np.append(linear_velocity_compose, angular_velocity_vector)
                    body.prescribed_velocity = velocity

                # Forward Euler method
                if self.numerical_method == "forward_euler":
                    location_new = body.location + linear_velocity_compose * dt
                    omega_axis_quaternion_dt = Quaternion.from_rotation(angular_velocity_vector * dt)
                    body.location = location_new
                    body.omega_axis_orientation = omega_axis_quaternion_dt * body.omega_axis_orientation
                    body.v0_axis = body.update_v0_axis_from_omega_axis(body.omega_axis_orientation)
                    velocity = np.append(linear_velocity_compose, angular_velocity_vector)
                    body.prescribed_velocity = velocity

                # Stochastic First Order
                if self.numerical_method == "stochastic_first_order":
                    random_rotation_vec_noise = np.random.randn(3)
                    stochastic_rotation_vec = np.sqrt(2 / self.gamma_r) * random_rotation_vec_noise * np.sqrt(dt)
                    omega_axis_quaternion_dt = Quaternion.from_rotation(angular_velocity_vector * dt +
                                                                        stochastic_rotation_vec)

                    body.omega_axis_orientation = omega_axis_quaternion_dt * body.omega_axis_orientation
                    omega_axis = body.update_omega_axis(body.omega_axis_orientation)
                    body.v0_axis = body.update_v0_axis_from_omega_axis(body.omega_axis_orientation)

                    intrinsic_swim_velocity = body.v0_axis * self.intrinsic_velocity[0]
                    random_translation = np.random.randn(3)
                    translational_noise_term = np.sqrt(2 / self.gamma_t) * random_translation * np.sqrt(dt)
                    linear_velocity_compose = intrinsic_swim_velocity + chem_prop
                    location_new = body.location + linear_velocity_compose * dt + translational_noise_term
                    body.location = location_new
                    velocity = np.append(linear_velocity_compose, angular_velocity_vector)
                    body.prescribed_velocity = velocity

                body.chem_surface_gradient = force_grad
                body.chem_torque_gradient = torque_grad

            else:
                # Use forward Euler method for the first step
                # chem_force = self.calc_tangential_grad_3D(self.body, *args, **kwargs)
                # chem_prop = self.mobility_alpha / (4 * np.pi) * chem_force

                if body.is_janus:
                    # Distribution case: returns force and torque gradients
                    force_grad, torque_grad = self.history_local_compose_3d_distribution(self.body, *args, **kwargs)
                    chem_prop = (self.mobility_alpha / (4 * np.pi)) * force_grad
                    chem_torque = (self.mobility_alpha / (4 * np.pi)) * torque_grad
                else:
                    # Point case: only returns a force gradient
                    force_grad = self.history_local_compose_3d_point(self.body, *args, **kwargs)
                    chem_prop = (self.mobility_alpha / (4 * np.pi)) * force_grad

                # Update position using Euler step
                intrinsic_swim_velocity = body.v0_axis * self.intrinsic_velocity[0]
                linear_velocity_compose = intrinsic_swim_velocity + chem_prop
                location_new = body.location + linear_velocity_compose * dt     # noise required
                body.location_new = location_new

                # Update orientation of omega axis and v_0 axis
                # noise required
                omega_axis = body.update_omega_axis(body.omega_axis_orientation)
                angular_velocity_vector = self.intrinsic_velocity[1] * omega_axis + chem_torque
                omega_axis_quaternion_dt = Quaternion.from_rotation(angular_velocity_vector * dt)
                body.omega_axis_orientation = omega_axis_quaternion_dt * body.omega_axis_orientation

                body.v0_axis = body.update_v0_axis_from_omega_axis(body.omega_axis_orientation)
                velocity = np.append(linear_velocity_compose, angular_velocity_vector)
                body.prescribed_velocity = velocity
                body.chem_surface_gradient = force_grad
                body.chem_torque_gradient = torque_grad

            # Update configuration
            body.location_history[step + 1, :] = location_new
            body.orientation = body.omega_axis_orientation.flip_self()
            body.orientation_history[step + 1, :] = body.omega_axis_orientation.flip_self()
            self.first_step = False
            self.velocities_previous_step = velocity

            return
