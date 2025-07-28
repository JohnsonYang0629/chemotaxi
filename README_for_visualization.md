# Chemical Concentration Field Visualization Tool

This script, `visualize_concentration.py`, is a post-processing tool designed to work with the output of the `chemotaxi` simulation package. It calculates and visualizes the chemical concentration field at a specific moment in time, based on the complete trajectory of a simulated particle.

## Features

-   Reads simulation parameters directly from the `*.inputfile`.
-   Loads particle trajectory data (`*.config` file).
-   Supports both **2D** and **3D** simulation domains.
    -   **For 2D**, it generates a 2D heatmap.
    -   **For 3D**, it now computes the concentration on a full 3D grid and generates a **3D scatter plot**, visualizing areas of significant concentration.
-   **Handles two types of 3D particles**:
    1.  **`non-janus`**: A simple point particle source.
    2.  **`janus`**: A spherical particle with a defined surface chemical distribution, accounting for both translation and rotation (via quaternions).
-   Uses `numba` for just-in-time parallel compilation to significantly accelerate calculations.
-   Optionally saves the calculated concentration grid to a file. For 3D, it saves to a more efficient binary `.npy` format.

## Prerequisites

-   Python 3.x
-   Numpy
-   Matplotlib (`mplot3d` toolkit included)
-   Numba
-   The `read_input.py` module from the original package.

Install libraries: `pip install numpy matplotlib numba`

## How to Use

### For 3D Simulations (Full 3D Visualization)

The primary update is for 3D visualization. The script no longer computes a 2D slice but a full 3D volume.

**Important Note on Performance:** Full 3D calculations are computationally intensive. The number of points to calculate is `resolution³`. A resolution of 50 (default) means 125,000 points. A resolution of 100 means 1,000,000 points. **Be cautious when increasing the `--resolution` for 3D simulations.**

**Example for a 3D `janus` simulation:**

1.  **Ensure your `*.inputfile` is configured for a Janus particle:**
    ```
    domain = 3D
    particle_type = janus
    surface_distribution_file = path/to/your/surface_sources.chem_dist.dat
    ```

2.  **Ensure your `*.config` trajectory file has 7 columns:** `x y z qw qx qy qz`.

3.  **Run the script:**
    ```bash
    python visualize_concentration.py \
        --input-file simulation_results/janus_3d.inputfile \
        --trajectory-file simulation_results/janus_3d.config \
        --time 100.0 \
        --output-file concentration_3d_janus.png \
        --resolution 60 \
        --save-dat
    ```
    -   This will create a 3D plot named `concentration_3d_janus.png`.
    -   It will also save the raw 3D concentration data to `concentration_3d_janus.npy`. You can load this file later using `data = np.load('concentration_3d_janus.npy')`.

## Command-Line Arguments

-   `--input-file` **(required)**: Path to the simulation input file (`*...inputfile`).
-   `--trajectory-file` **(required)**: Path to the trajectory file (`*...config`).
-   `--time` **(required)**: The simulation time `t1` for the snapshot.
-   `--output-file` (optional): Filename for the saved plot. Default: `concentration_field.png`.
-   `--resolution` (optional): The number of points per axis for the grid. **Use with caution for 3D.** Default: `50`.
-   `--padding` (optional): Space to add around the trajectory bounds. Default: `5.0`.
-   `--save-dat` (optional flag): If specified, saves the concentration grid. For 3D, this creates a `.npy` file.