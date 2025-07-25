import argparse
import sys

import numpy as np
from shutil import copyfile
from functools import partial
import time
import pickle

# Find project functions
from read_input import read_input
from read_input import read_vertex_file
from read_input import read_chem_dist_file
from body import body_2D
from body import body_3D
from integrator.integrator_2D import ChemoIntegrator2D
from integrator.integrator_3D import ChemoIntegrator3D
import chem_functions


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # Get command line arguments
    parser = argparse.ArgumentParser(description='Run a chemotaxi simulation and save trajectory.')
    parser.add_argument('--input-file', dest='input_file', type=str, default='data.main', help='name of the input file')
    args = parser.parse_args()
    input_file = args.input_file

    # Read input file
    read = read_input.ReadInput(input_file)

    scheme = read.scheme
    domain = read.domain
    numerical_method = read.numerical_method
    n_steps = read.n_steps
    n_save = read.n_save
    dt = read.dt

    intrinsic_linear_velocity = read.intrinsic_linear_velocity
    intrinsic_angular_velocity = read.intrinsic_angular_velocity
    intrinsic_velocity = np.array([intrinsic_linear_velocity, intrinsic_angular_velocity])

    output_name = read.output_name
    copyfile(input_file, output_name + '.inputfile')

    structure_file_name = read.structure
    structure_ref_config = read_vertex_file.read_vertex_file(structure_file_name[0])

    if read.random_state != 'None':
        random_state_file = read.random_state
        print(f"Attempting to load random state from: {random_state_file}")
        try:
            with open(random_state_file, 'rb') as f:
                # load random state from file
                state = pickle.load(f)
                # apply random state to numpy
                np.random.set_state(state)
            print(f"Successfully loaded random state from '{random_state_file}'.")
        except FileNotFoundError:
            print(f"ERROR: The specified random_state file '{random_state_file}' was not found.")
            exit(1)
        except Exception as e:
            print(f"ERROR: An error occurred while loading the random state: {e}")
            exit(1)

    if read.numerical_method == "stochastic_first_order":
        output_state_filename = f"{read.output_name}.random_state"
        print(f"Saving the initial random state for this run to '{output_state_filename}'.")
        try:
            # get and write current numpy random state
            current_state = np.random.get_state()
            with open(output_state_filename, 'wb') as f:
                pickle.dump(current_state, f)
        except IOError as e:
            print(f"ERROR: Could not save the random state to '{output_state_filename}': {e}")

    # Create droplet body
    if domain == '2D':
        initial_struct_location_2D = read.initial_position_2D
        initial_struct_orientations_2D = read.initial_orientation_2D_vector
        body = body_2D.Body2D(initial_struct_location_2D, initial_struct_orientations_2D, n_steps)
        body.location_history[0, :] = initial_struct_location_2D
        integrator = ChemoIntegrator2D(body, scheme, domain, numerical_method)
    elif domain == '3D':
        initial_struct_location_3D = read.initial_position_3D
        initial_omega_axis_orientations_3D = read.initial_orientation_3D_quaternion
        body = body_3D.Body3D(initial_struct_location_3D, initial_omega_axis_orientations_3D, structure_ref_config, n_steps)
        body.location_history[0, :] = initial_struct_location_3D
        body.orientation_history[0, :] = initial_omega_axis_orientations_3D.flip_self()
        body.orientation = initial_omega_axis_orientations_3D.flip_self()
        integrator = ChemoIntegrator3D(body, scheme, domain, numerical_method)
        chemical_dist_file = read.chemical_distribution_file
        sigma_values = read_chem_dist_file.read_chemical_distribution_file(chemical_dist_file[0])
        body.set_sigma_distribution(sigma_values)

    else:
        print('Domain should use \"2D\" or \"3D\". \n')
        exit()

    integrator.peclet_number = read.peclet_number
    integrator.mobility_alpha = read.mobility_alpha
    integrator.intrinsic_velocity = intrinsic_velocity
    integrator.gamma_r = read.gamma_r
    integrator.gamma_t = read.gamma_t

    if domain == '2D':
        integrator.rotation_matrix_2d = chem_functions.rotation_matrix_2d
        integrator.calc_surface_gradient_circle = partial(chem_functions.calc_surface_gradient_circle,
                                                          acceleration=read.acceleration,
                                                          core=read.core,
                                                          peclet_number=read.peclet_number,
                                                          structure_ref_config=structure_ref_config,
                                                          dt=dt)
    elif domain == '3D':
        if read.particle_type == 'non_janus':
            integrator.history_local_compose_3d_point = partial(chem_functions.history_local_compose_3d_point,
                                                                acceleration=read.acceleration,
                                                                core=read.core,
                                                                peclet_number=read.peclet_number,
                                                                structure_ref_config=structure_ref_config,
                                                                dt=dt)
        elif read.particle_type == 'janus':
            integrator.history_local_compose_3d_distribution = partial(chem_functions.history_local_compose_3d_distribution,
                                                                       acceleration=read.acceleration,
                                                                       core=read.core,
                                                                       peclet_number=read.peclet_number,
                                                                       structure_ref_config=structure_ref_config,
                                                                       dt=dt)
        elif read.particle_type == 'default':
            integrator.calc_tangential_grad_3D = partial(chem_functions.calc_tangential_grad_3D,
                                                         acceleration=read.acceleration,
                                                         core=read.core,
                                                         peclet_number=read.peclet_number,
                                                         structure_ref_config=structure_ref_config,
                                                         dt=dt)

    # Loop over time steps
    start_time = time.time()
    if read.save_clones == 'one_file':
        buffering = max(1,  n_steps // n_save // 200)
        output_file_name = output_name + '.config'
        loc_file = open(output_file_name, 'w', buffering=buffering)
        velocity_file_name = output_name + '.velocity.dat'
        velocity_file = open(velocity_file_name, 'w', buffering=buffering)
        chem_force_file_name = output_name + '.chemforce.dat'
        chem_force_file = open(chem_force_file_name, 'w', buffering=buffering)
        time_log_file_name = output_name + '.time.log'
        time_log_file = open(time_log_file_name, 'w', buffering=buffering)

    for step in range(read.initial_step, n_steps):
        # Save data if...
        if (step % n_save) == 0 and step >= 0:
            elapsed_time = time.time() - start_time
            print('Step = ', step, ', wallclock time = ', elapsed_time)

            if domain == '2D':
                loc_file.write('%s %s\n' % (body.location[0], body.location[1]))
                velocity_file.write('%s %s %s\n' % (body.prescribed_velocity[0],
                                                    body.prescribed_velocity[1],
                                                    body.prescribed_velocity[2]))
                chem_force_file.write('%s %s\n' % (body.chem_surface_gradient[0], body.chem_surface_gradient[1]))
                time_log_file.write(str(elapsed_time) + '\n')
            elif domain == '3D':
                loc_file.write('%s %s %s %s %s %s %s\n' % (body.location[0],
                                                           body.location[1],
                                                           body.location[2],
                                                           body.orientation[3],
                                                           body.orientation[0],
                                                           body.orientation[1],
                                                           body.orientation[2]))
                velocity_file.write('%s %s %s %s %s %s\n' % (body.prescribed_velocity[0],
                                                             body.prescribed_velocity[1],
                                                             body.prescribed_velocity[2],
                                                             body.prescribed_velocity[3],
                                                             body.prescribed_velocity[4],
                                                             body.prescribed_velocity[5]))
                chem_force_file.write('%s %s %s %s %s %s\n' % (body.chem_surface_gradient[0],
                                                               body.chem_surface_gradient[1],
                                                               body.chem_surface_gradient[2],
                                                               body.chem_torque_gradient[0],
                                                               body.chem_torque_gradient[1],
                                                               body.chem_torque_gradient[2]))
                time_log_file.write(str(elapsed_time) + '\n')

        integrator.advance_time_step(dt, step=step)

    # Save final data if...
    if ((step + 1) % n_save) == 0 and step >= 0:
        elapsed_time = time.time() - start_time
        print('Step = ', step + 1, ', wallclock time = ', elapsed_time)
        if domain == '2D':
            loc_file.write('%s %s\n' % (body.location[0], body.location[1]))
            velocity_file.write('%s %s %s\n' % (body.prescribed_velocity[0],
                                                body.prescribed_velocity[1],
                                                body.prescribed_velocity[2]))
            chem_force_file.write('%s %s\n' % (body.chem_surface_gradient[0], body.chem_surface_gradient[1]))
            time_log_file.write(str(elapsed_time) + '\n')
        elif domain == '3D':
            loc_file.write('%s %s %s %s %s %s %s\n' % (body.location[0],
                                                       body.location[1],
                                                       body.location[2],
                                                       body.orientation[3],
                                                       body.orientation[0],
                                                       body.orientation[1],
                                                       body.orientation[2]))
            velocity_file.write('%s %s %s %s %s %s\n' % (body.prescribed_velocity[0],
                                                         body.prescribed_velocity[1],
                                                         body.prescribed_velocity[2],
                                                         body.prescribed_velocity[3],
                                                         body.prescribed_velocity[4],
                                                         body.prescribed_velocity[5]))
            chem_force_file.write('%s %s %s %s %s %s\n' % (body.chem_surface_gradient[0],
                                                           body.chem_surface_gradient[1],
                                                           body.chem_surface_gradient[2],
                                                           body.chem_torque_gradient[0],
                                                           body.chem_torque_gradient[1],
                                                           body.chem_torque_gradient[2]))
            time_log_file.write(str(elapsed_time) + '\n')

    with open(output_name + '.time', 'w') as f:
        f.write(str(time.time() - start_time) + '\n')
