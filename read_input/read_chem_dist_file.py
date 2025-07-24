'''
Small module to read a chemical substance distribution file of the rigid bodies.
'''
import numpy as np


def read_chemical_distribution_file(name_file):
    """
    Reads a chemical distribution file.
    The file should contain one column of numbers, where each row corresponds
    to a node in the vertex file.

    Args:
        name_file (str): The path to the .chem_dist.dat file.

    Returns:
        np.ndarray or None: A numpy array with shape (N_nodes, 1) containing
                             the concentration values, or None if the file
                             name is 'None'.
    """
    # If the file name is not provided in the input file, do nothing.
    if name_file is None or name_file == 'None':
        return None

    try:
        # Use numpy.loadtxt for efficient reading of numerical data.
        concentrations = np.loadtxt(name_file)

        # Ensure the output is a column vector (N, 1) for consistent matrix operations later.
        return concentrations.reshape(-1, 1)

    except IOError:
        print(f"Error: Could not read the chemical distribution file at '{name_file}'. Please check the path.")
        return None
    except ValueError:
        print(f"Error: The file '{name_file}' contains non-numerical data. Please check the file content.")
        return None