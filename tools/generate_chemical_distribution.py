import numpy as np
import os


def generate_chemical_distribution(input_path, input_filename):
    """
    Reads sphere vertex data, calculates a chemical distribution based on
    hemisphere, and writes the distribution to a new file.

    This function mimics the logic of the original MATLAB script.
    - A value of 1 is assigned to points in the northern hemisphere (z >= 0).
    - A value of 0 is assigned to points in the southern hemisphere (z < 0).

    Args:
        input_path (str): The directory path where the input file is located.
        input_filename (str): The name of the input vertex file.
    """
    full_input_path = os.path.join(input_path, input_filename)

    # --- 2. Read the data from the file ---
    try:
        print(f"Reading data from: {full_input_path}")
        # np.loadtxt is a robust way to read space-delimited data.
        # We don't need to read the number of points from the first line
        # as numpy can determine the size of the array automatically.
        # If the file truly has a header line with the count, use skiprows=1.
        # Based on the MATLAB script, it seems the first line is NOT part of the coordinates.
        with open(full_input_path, 'r') as f:
            # The MATLAB script reads the number of points, but it's not strictly
            # necessary if the file only contains coordinate data. For a direct
            # translation, we can read it but won't use it.
            # If your file format is simply a list of coordinates, you can use:
            # coordinates = np.loadtxt(full_input_path)

            # Assuming the first line is a header to be skipped.
            next(f, None)  # Skips the first line (number of points)
            coordinates = np.loadtxt(f)

        num_points = coordinates.shape[0]
        print(f"Successfully read coordinates for {num_points} points.")

    except FileNotFoundError:
        print(f"Error: The file was not found at {full_input_path}")
        return
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return

    # --- 3. Calculate the chemical distribution ---
    # In spherical coordinates, the polar angle theta is measured from the positive Z-axis.
    # theta in [0, pi/2] corresponds to the northern hemisphere (z >= 0).
    # theta in (pi/2, pi] corresponds to the southern hemisphere (z < 0).
    # A simpler and more direct check is to just use the sign of the Z coordinate.

    # Extract the Z coordinates (the 3rd column, index 2)
    z_coords = coordinates[:, 2]

    # Initialize the distribution array with zeros
    chem_distribution = np.zeros(num_points, dtype=int)

    # Set points in the northern hemisphere (where z >= 0) to 1
    # This is much more efficient than calculating acos.
    is_northern_hemisphere = z_coords >= 0
    chem_distribution[is_northern_hemisphere] = 1

    print("Chemical distribution calculated successfully.")

    # --- 4. Write the distribution to a .chem_dist.dat file ---
    # Construct the output filename
    name, _ = os.path.splitext(input_filename)
    output_filename = name + '.chem_dist.dat'
    full_output_path = os.path.join(input_path, output_filename)

    try:
        # Save the array to a text file.
        # fmt='%d' ensures the output is integer format.
        np.savetxt(full_output_path, chem_distribution, fmt='%d')
        print(f"Successfully generated file: {full_output_path}")

    except Exception as e:
        print(f"An error occurred while writing the file: {e}")
        return

    print("\nScript finished successfully!")


if __name__ == '__main__':
    # --- Configuration ---
    # You can change these paths to match your file locations.
    # For better practice, consider using command-line arguments (e.g., with argparse).
    # NOTE: The path is an example. Please adjust it to your actual directory.
    path_name = "../chemotaxi/structures"
    file_name = "sphere_N_42_R_1.vertex"

    # Run the main function
    generate_chemical_distribution(path_name, file_name)
