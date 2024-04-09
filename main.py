import argparse
import numpy as np
from shutil import copyfile
import time

# Find project functions
from read_input import read_input
from body import body_2D
from body import body_3D
from integrator.integrator_2D import ChemoIntegrator2D
from integrator.integrator_3D import ChemoIntegrator3D

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
    n_steps = read.n_steps
    n_save = read.n_save
    dt = read.dt

    initial_struct_location = read.initial_position
    initial_struct_orientations = read.initial_orientation
    intrinsic_linear_velocity = read.intrinsic_linear_velocity
    intrinsic_angular_velocity = read.intrinsic_angular_velocity
    intrinsic_velocity = np.array([intrinsic_linear_velocity, intrinsic_angular_velocity])

    output_name = read.output_name
    copyfile(input_file, output_name + '.inputfile')

    # Create droplet body
    if domain == '2D':
        body = body_2D.Body_2D(initial_struct_location, initial_struct_orientations)
        integrator = ChemoIntegrator2D(body, scheme, domain)
    elif domain == '3D':
        body = body_3D.Body_3D(initial_struct_location, initial_struct_orientations)
        integrator = ChemoIntegrator3D(body, scheme, domain)
    else:
        print('Domain should use \"2D\" or \"3D\". \n')
        exit()



    # Loop over time steps
    start_time = time.time()
    if read.save_clones == 'one_file':
        buffering = max(1,  n_steps // n_save // 200)
        output_file_name = output_name + '.config'
        open(output_file_name, 'w', buffering=buffering)
        velocity_file_name = output_name + '.velocity.dat'
        open(velocity_file_name, 'w', buffering=buffering)
        chem_force_file_name = output_name + '.chemforce.dat'
        open(chem_force_file_name, 'w', buffering=buffering)

    for step in range(read.initial_step, n_steps):
        # Save data if...
        if (step % n_save) == 0 and step >= 0:
            elapsed_time = time.time() - start_time
            print('Step = ', step, ', wallclock time = ', time.time() - start_time)


        integrator.advance_time_step(dt, step=step)

        # Save final data if...
        if ((step + 1) % n_save) == 0 and step >= 0:
