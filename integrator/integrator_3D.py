import numpy as np
from integrator.quaternion import Quaternion

class ChemoIntegrator3D(object):

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
        self.intrinsic_velocity = np.array([0, 0])  # compact vector [v_0,omega_0]
        self.gamma_r = 0.0
        self.gamma_t = 0.0

        # Optional variables
        self.calc_tangential_grad_3D = None
        self.rotation_matrix_3d = None

    def advance_time_step(self, dt, *args, **kwargs):
        """
        Advance time step with integrator self.scheme
        """
        return getattr(self, self.scheme)(dt, *args, **kwargs)

    def noise(self, dt, *args, **kwargs):
        return

    def history_local_compose_3d(self, dt, *args, **kwargs):
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
                chem_force = self.calc_tangential_grad_3D(self.body, *args, **kwargs)
                chem_prop = self.mobility_alpha / (4 * np.pi) * chem_force

                intrinsic_swim_velocity = body.v0_axis * self.intrinsic_velocity[0]
                linear_velocity_compose = intrinsic_swim_velocity + chem_prop

                omega_axis = body.update_omega_axis(body.omega_axis_orientation)
                angular_velocity_vector = self.intrinsic_velocity[1] * omega_axis

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

                body.chem_surface_gradient = chem_force

            else:
                # Use forward Euler method for the first step
                chem_force = self.calc_tangential_grad_3D(self.body, *args, **kwargs)
                chem_prop = self.mobility_alpha / (4 * np.pi) * chem_force

                # Update position using Euler step
                intrinsic_swim_velocity = body.v0_axis * self.intrinsic_velocity[0]
                linear_velocity_compose = intrinsic_swim_velocity + chem_prop
                location_new = body.location + linear_velocity_compose * dt     # noise required
                body.location_new = location_new

                # Update orientation of omega axis and v_0 axis
                # noise required
                omega_axis = body.update_omega_axis(body.omega_axis_orientation)
                angular_velocity_vector = self.intrinsic_velocity[1] * omega_axis
                omega_axis_quaternion_dt = Quaternion.from_rotation(angular_velocity_vector * dt)
                body.omega_axis_orientation = omega_axis_quaternion_dt * body.omega_axis_orientation

                body.v0_axis = body.update_v0_axis_from_omega_axis(body.omega_axis_orientation)
                velocity = np.append(linear_velocity_compose, angular_velocity_vector)
                body.prescribed_velocity = velocity
                body.chem_surface_gradient = chem_force

            # Update configuration
            body.location_history[step + 1, :] = location_new
            self.first_step = False
            self.velocities_previous_step = velocity

            return
