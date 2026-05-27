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
import chem_solver


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

    # structure_file_name = read.structure
    # structure_ref_config = read_vertex_file.read_vertex_file(structure_file_name[0])

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

    # Create a list to hold all body objects
    bodies = []

    if domain == '2D':
        for i in range(read.droplet_num):
            location = read.initial_positions_2D[i]
            orientation = read.initial_orientations_2D[i]
            structure_file = read.structures[i][0] if i < len(read.structures) else read.structures[0][0]
            structure_ref_config = read_vertex_file.read_vertex_file(structure_file)
            body = body_2D.Body2D(location, orientation, structure_ref_config, n_steps)
            body.location_history[0, :] = location
            body.ID = i
            body.peclet_number = read.peclet_numbers[i]
            body.mobility_alpha = read.mobility_alphas[i]

            particle_type = read.particle_types[i] if i < len(read.particle_types) else 'non_janus'
            if particle_type == 'janus':
                chem_dist_file = read.chemical_distribution_files[i] \
                    if i < len(read.chemical_distribution_files) else 'None'
                if chem_dist_file != 'None':
                    sigma_values = read_chem_dist_file.read_chemical_distribution_file(chem_dist_file)
                    body.set_sigma_distribution(sigma_values)
                else:
                    body.set_sigma_distribution(None)
            else:
                body.is_janus = False

            mob_dist_file = read.mobility_distribution_files[i] if i < len(read.mobility_distribution_files) else 'None'
            if mob_dist_file != 'None':
                mob_values = read_chem_dist_file.read_chemical_distribution_file(mob_dist_file)
                body.set_mobility_distribution(mob_values)
            else:
                body.set_mobility_distribution(None)

            bodies.append(body)

    elif domain == '3D':
        for i in range(read.droplet_num):
            location = read.initial_positions_3D[i]
            orientation_quat = read.initial_orientations_3D_quaternion[i]
            structure_file = read.structures[i][0] if i < len(read.structures) else read.structures[0][0]
            structure_ref_config = read_vertex_file.read_vertex_file(structure_file)
            body = body_3D.Body3D(location, orientation_quat, structure_ref_config, n_steps)
            body.location_history[0, :] = location
            body.orientation_history[0, :] = orientation_quat.flip_self()
            body.orientation = orientation_quat.flip_self()
            body.ID = i
            body.peclet_number = read.peclet_numbers[i]
            body.mobility_alpha = read.mobility_alphas[i]

            particle_type = read.particle_types[i] if i < len(read.particle_types) else 'non_janus'
            if particle_type == 'janus':
                chem_dist_file = read.chemical_distribution_files[i] \
                    if i < len(read.chemical_distribution_files) else 'None'
                if chem_dist_file != 'None':
                    sigma_values = read_chem_dist_file.read_chemical_distribution_file(chem_dist_file)
                    body.set_sigma_distribution(sigma_values)
                else:
                    body.set_sigma_distribution(None)
            else:
                body.is_janus = False

            bodies.append(body)

    else:
        print('Domain should use \"2D\" or \"3D\". \n')
        exit()

    # --- INTEGRATOR INITIALIZATION ---
    # Initialize the integrator by passing the entire list of bodies.
    # The integrator will now manage the state of these bodies.
    if domain == '2D':
        integrator = ChemoIntegrator2D(bodies, scheme, domain, numerical_method)
    elif domain == '3D':
        integrator = ChemoIntegrator3D(bodies, scheme, domain, numerical_method)

    integrator.intrinsic_velocity = intrinsic_velocity
    integrator.gamma_r = read.gamma_r
    integrator.gamma_t = read.gamma_t

    if domain == '2D':
        integrator.rotation_matrix_2d = chem_solver.rotation_matrix_2d
        integrator.history_local_compose_2d_multi_body = partial(chem_solver.history_local_compose_2d_multi_body,
                                                                 dt=dt)
    elif domain == '3D':
        integrator.history_local_compose_3d_multi_body = partial(chem_solver.history_local_compose_3d_multi_body,
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

            if step % 100 == 0:
                time_log_file.write(f'{step} {elapsed_time}\n')

            loc_file.write(str(read.droplet_num) + '\n')
            velocity_file.write(str(read.droplet_num) + '\n')
            chem_force_file.write(str(read.droplet_num) + '\n')

            if domain == '2D':
                for body in integrator.bodies:
                    loc_file.write('%s %s %s %s\n' % (body.location[0],
                                                      body.location[1],
                                                      body.orientation[0],
                                                      body.orientation[1]))
                    velocity_file.write('%s %s %s\n' % (body.prescribed_velocity[0],
                                                        body.prescribed_velocity[1],
                                                        body.prescribed_velocity[2]))

                    chem_force_file.write('%s %s %s\n' % (body.chem_surface_gradient[0],
                                                          body.chem_surface_gradient[1],
                                                          body.chem_torque_gradient))

            elif domain == '3D':
                for body in integrator.bodies:
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

        integrator.advance_time_step(dt, step=step)

    # Save final data if...
    if ((step + 1) % n_save) == 0 and step >= 0:
        elapsed_time = time.time() - start_time
        print('Step = ', step + 1, ', wallclock time = ', elapsed_time)
        time_log_file.write(f'{step} {elapsed_time}\n')

        loc_file.write(str(read.droplet_num) + '\n')
        velocity_file.write(str(read.droplet_num) + '\n')
        chem_force_file.write(str(read.droplet_num) + '\n')

        if domain == '2D':
            for body in integrator.bodies:
                loc_file.write('%s %s %s %s\n' % (body.location[0],
                                                  body.location[1],
                                                  body.orientation[0],
                                                  body.orientation[1]))
                velocity_file.write('%s %s %s\n' % (body.prescribed_velocity[0],
                                                    body.prescribed_velocity[1],
                                                    body.prescribed_velocity[2]))
                chem_force_file.write('%s %s %s\n' % (body.chem_surface_gradient[0],
                                                      body.chem_surface_gradient[1],
                                                      body.chem_torque_gradient))

        elif domain == '3D':
            for body in integrator.bodies:
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

    with open(output_name + '.time', 'w') as f:
        f.write(str(time.time() - start_time) + '\n')
