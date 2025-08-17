# Chemical Concentration Field Visualization Tool

This script, `visualize_concentration.py`, is a post-processing tool designed to work with the output of the chemotaxi simulation package. 
It calculates and visualizes the chemical concentration field at a specific moment in time, based on the complete trajectory of the simulated particles.

## Features

- Reads simulation parameters directly from the `*.inputfile`, 
including support for individual parameters (like Peclet numbers) for multiple particles.
- Loads and parses the new multi-particle trajectory data (`*.config` file).
- Supports both **2D** and **3D** simulation domains.
    - **2D Mode**: Generates static heatmaps or full animations of the concentration evolution.
    - **3D Mode**: Efficiently computes and renders two key 2D concentration slices, 
visualizing them in context with the 3D trajectories, avoiding expensive calculations of the full 3D volume.
- **Handles two types of 3D particles**:
    1.  **`non-janus`**: A simple point particle source.
    2.  **`janus`**: A spherical particle with a defined surface chemical distribution, accounting for both translation and rotation (via quaternions).
- Uses `numba` for just-in-time parallel compilation to significantly accelerate calculations.
- Flexible Command-Line Interface:
  1. Supports using a file prefix for simplified input, automatically matching `.inputfile` and `.config` files.
  2. Provides multiple 3D slice positioning modes, including auto-centering, targeting the end of a trajectory, or specifying exact coordinates manually.
- Optionally saves the calculated concentration grid to a file. For 3D, it saves to a more efficient binary `.npy` format.

## Prerequisites

-   Python 3.x
-   Numpy
-   Matplotlib (`mplot3d` toolkit included)
-   Numba
-   The `read_input.py` module from the original package.

Install libraries: `pip install numpy matplotlib numba`

## How to Use

### For 2D Simulations (Static and Animated Outputs)
#### Example 1: Generate a 2D Static Plot
To create a static image of the concentration field for a 2D simulation at time `t=100.0`.
```
python visualize_concentration.py \
    --input-file simulation_results/chemo_2d_test.inputfile \
    --trajectory-file simulation_results/chemo_2d_test.config \
    --mode static \
    --time 100.0 \
    --output-file static_2d_concentration
```
This will generate `static_2d_concentration.png`.

#### Example 2: Generate a 2D Animation (GIF format)
To create a full animation for a 2D simulation, with each frame corresponding to `1.0` units of simulation time.
```
python visualize_concentration.py \
    --input-file simulation_results/chemo_2d_test.inputfile \
    --trajectory-file simulation_results/chemo_2d_test.config \
    --mode animation \
    --frame-interval 1.0 \
    --writer gif \
    --output-file animation_2d_chemo
```
This will generate `animation_2d_chemo.gif`.

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

- `--input-file` **(required)**: Path to the simulation input file (`*...inputfile`).
- `--trajectory-file` **(required)**: Path to the trajectory file (`*...config`).
- `--output-file` (optional): Filename for the saved plot (without extension). Default: `concentration_field.png`.
- `--resolution` (optional): The number of points per axis for the grid. **Use with caution for 3D.** Default: `50`.
- `--mode` (optional):The operating mode: `static` or `animation`. Default: `static`.
- `--time` **(required for static mode)**: The simulation time `t1` for the snapshot. 
- `frame-interval` (optional for animation): The amount of simulation time between each frame of the animation. Default: `5.0`.
- `writer` (optional for animation): The writer to use for saving animations: gif or ffmpeg. Default: `gif`.
- `--padding` (optional): Space to add around the trajectory bounds. Default: `5.0`.
- `--save-dat` (optional for static): If specified, saves the concentration grid. For 3D, this creates a `.npy` file.